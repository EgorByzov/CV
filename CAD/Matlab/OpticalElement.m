classdef (Abstract) OpticalElement < IRhinoExportable & ISTLExportable & ICloneable & IDrawable & ISerializable
    
    properties (Access = protected)
        surfacesListeners = {};
    end
    
    properties (Access = {?OpticalElement, ?Serializator})        
        surfaces;
        cs;
    end % properties
    
    properties (Dependent)
        Surfaces;
        NumSurfaces;
        CS;
        IsExtrudeOn;
    end % properties
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        function value = get.CS(obj)
            value = obj.cs;
        end
        function set.CS(obj, cs)
            ValueChecker.checkClass(cs, 'CoordinateSystem');
            ValueChecker.checkAndUpdateSize(cs, [1 1]);
            
            if ~isequal(cs, obj.cs)
                for i = 1:obj.NumSurfaces
                    curCS = obj.surfaces{i}.CS;
                    
                    % get surface CS in oldCS
                    center = GeometryConvertor.getLocalPoint(curCS.Center, obj.cs);
                    axY = GeometryConvertor.getLocalVector(curCS.AxisY, obj.cs);
                    axZ = GeometryConvertor.getLocalVector(curCS.AxisZ, obj.cs);
                    
                    % get surface CS in global CS from new CS
                    curCS.Center = GeometryConvertor.getGlobalPoint(center, cs);
                    curCS.setAxes( GeometryConvertor.getGlobalVector(axY, cs),...
                        GeometryConvertor.getGlobalVector(axZ, cs) );
                end
                
                obj.cs = cs;
            end
        end % function
        
        function value = get.Surfaces(obj)
            value = obj.surfaces;
        end % function
        function set.Surfaces(obj, surfaces)
            ValueChecker.checkClass(surfaces, 'Surface');
            obj.setSurfaces(surfaces);
        end % function
        
        function value = get.NumSurfaces(obj)
            value = length(obj.surfaces);
        end % function
        
        function value = get.IsExtrudeOn(obj)
            if (obj.isAxisymmetrical)
                value = obj.surfaces{1}.IsProfileExtruded;
            else
                value = false;
            end            
        end % function
        function set.IsExtrudeOn(obj, value)
            if (obj.isAxisymmetrical && obj.IsExtrudeOn ~= value)
                for i = 1 : obj.NumSurfaces
                    obj.surfaces{i}.IsProfileExtruded = value;
                end % for
            end
        end % function
        
    end % methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Public %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access = public)
        
        function res = isAxisymmetrical(obj)
            res = true;
            if isa(obj, 'LensCurveExtrude')
                return;
            end
            for i = 1 : obj.NumSurfaces
                if ~isa(obj.surfaces{i}, 'AxisymSurface')
                    res = false;
                    return;
                end % if
            end % for
        end % function
        
        function res = isOpticalElementOptimizationAvailable(obj)
            res = false;
            
            if isa(obj, 'IOptimizable')
                res = obj.isOptimizationAvailable();
            else
                for i = 1:length(obj.surfaces)
                    if isa(obj.surfaces{i}, 'IOptimizable') && obj.surfaces{i}.isOptimizationAvailable()
                        res = true;
                        return;
                    end
                end
            end
        end % function
        
        function res = isOpticalElementOptimizable(obj)
            res = false;
            
            if isa(obj, 'IOptimizable')
                res = true;
            else
                for i = 1:length(obj.surfaces)
                    if isa(obj.surfaces{i}, 'IOptimizable')
                        res = true;
                        return;
                    end
                end
            end
        end % function
        
        function res = getRhinoCommandsString(obj, version, params)
            res = [];
            for i = 1 : obj.NumSurfaces
                curStr = obj.surfaces{i}.getRhinoCommandsString(version, params);
                res = [res curStr];
            end % for
            
            if isa(obj, 'SolidElement') && ~isempty(res)
                res = [res  '_SelAll\n_Join\n'];
            end
        end % function
        
        function res = getRhinoCommandsStringMoldable(obj, params)
            res = obj.getRhinoCommandsStringMoldable(params);
        end % function
        
        function surfs = getPrimitiveSets(obj, approximOptions)
            surfs = cell(1, length(obj.surfaces));
            for i = 1:length(obj.surfaces)
                surfs{i} = obj.surfaces{i}.getPrimitiveSet(approximOptions);
            end % for
        end % function
        
        function rawTriangs = getTriangles(obj, apprOpt)
            surfs = obj.surfaces;
            numSurfs = obj.NumSurfaces;
            triangSets = cell(numSurfs, 1);
            rawTriangs = [];
            
            % generate rays
            tmpRays = zeros(numSurfs, 7);
            for i = 1 : numSurfs
                % create triangle sets in such manner to have real normals
                % of triangles (not computed using special function)
                triangSets{i} = TriangleSet(surfs{i}.getTriangles(apprOpt), surfs{i}.CS);
                tmpRays(i, :) = triangSets{i}.getOutgoingRay();
            end % for
            
            % compute number of intersections with every solid's surface
            numInts = zeros(numSurfs, numSurfs);
            for i = 1 : numSurfs
                numInts(:, i) = triangSets{i}.getNumInters(tmpRays);
            end % for
            numInts = sum(numInts, 2);
            
            for i = 1 : numSurfs
                curTriangles = triangSets{i}.Triangles;
                if mod(numInts(i), 2) == 0
                    tmp = curTriangles(:, 4:6);
                    curTriangles(:, 4:6) = curTriangles(:, 7:9);
                    curTriangles(:, 7:9) = tmp;
                    curTriangles(:, 10:12) = -curTriangles(:, 10:12);
                end
                
                % convert coordinates to optical element`s CS
                vectors = GeometryConvertor.getGlobalVector(curTriangles(:, 10:12), triangSets{i}.CS);
                curTriangles(:, 10:12) = GeometryConvertor.getLocalVector(vectors, obj.cs);
                
                x = [curTriangles(:, 1) curTriangles(:, 4) curTriangles(:, 7)];
                y = [curTriangles(:, 2) curTriangles(:, 5) curTriangles(:, 8)];
                z = [curTriangles(:, 3) curTriangles(:, 6) curTriangles(:, 9)];
                defSize = size(x);
                
                x = reshape(x, [numel(x), 1]);
                y = reshape(y, [numel(y), 1]);
                z = reshape(z, [numel(z), 1]);
                
                points = GeometryConvertor.getGlobalPoint([x y z], triangSets{i}.CS);
                points = GeometryConvertor.getLocalPoint(points, obj.cs);
                
                x = reshape(points(:,1), defSize);
                y = reshape(points(:,2), defSize);
                z = reshape(points(:,3), defSize);
                
                curTriangles(:,1) = x(:,1); curTriangles(:,4) = x(:,2); curTriangles(:,7) = x(:,3);
                curTriangles(:,2) = y(:,1); curTriangles(:,5) = y(:,2); curTriangles(:,8) = y(:,3);
                curTriangles(:,3) = z(:,1); curTriangles(:,6) = z(:,2); curTriangles(:,9) = z(:,3);
                
                rawTriangs = [rawTriangs; curTriangles];
            end % for
        end % function
        
        function [pointsCloud, lines] = getPointsCloudAndLines(obj)
            pointsCloud = cell(1, obj.NumSurfaces);
            lines = cell(1, obj.NumSurfaces);
            
            for j = 1 : obj.NumSurfaces
                [pointsCloud{j}, lines{j}] = obj.surfaces{j}.getPointsCloudAndLines();
            end % for
        end % function
        
        function line2D = getProfile2D(obj, phi, numPoints)
            if nargin < 2
                phi = 0;
            end % if
            if nargin < 3
                numPoints = 200;
            end % if
            
            % construct profile and pack
            numSrfs = obj.NumSurfaces;
            lines = cell(numSrfs, 1);
            for i = 1 : numSrfs
                lines{i} = obj.surfaces{i}.getProfile2D(phi, numPoints);
            end % for
            line2D = MultiLine2D(lines);
        end % function
        
        function line2D = getSection2D(obj, phi, numPoints)
            if nargin < 2
                phi = 0;
            end % if
            if nargin < 3
                numPoints = 200;
            end % if
            
            % construct profile and pack
            numSrfs = obj.NumSurfaces;
            lines = cell(numSrfs, 1);
            for i = 1 : numSrfs
                lines{i} = obj.surfaces{i}.getSection2D(phi, numPoints);
            end % for
            line2D = MultiLine2D(lines);
        end % function
        
        function box = getElementBox(obj)
            pointsCloud = obj.getPointsCloudAndLines();
            xMin = inf; yMin = inf; zMin = inf;
            xMax = -inf; yMax = -inf; zMax = -inf;
            
            for i = 1:length(pointsCloud)
                curPointsCloud = pointsCloud{i};
                xMin = min( min(min( curPointsCloud(:,:,1) )), xMin );
                xMax = max( max(max( curPointsCloud(:,:,1) )), xMax );
                
                yMin = min( min(min( curPointsCloud(:,:,2) )), yMin );
                yMax = max( max(max( curPointsCloud(:,:,2) )), yMax );
                
                zMin = min( min(min( curPointsCloud(:,:,3) )), zMin );
                zMax = max( max(max( curPointsCloud(:,:,3) )), zMax );
            end
            box = [xMin xMax;
                yMin yMax;
                zMin zMax];
        end
        
        function res = isRegistratorOutside(obj, registrator)
            ValueChecker.checkClass(registrator, 'Surface');
            ValueChecker.checkAndUpdateSize(registrator, [1 1]);
            ValueChecker.checkConstraint(registrator.IsRegistrator, '==', true);
            
            % get element`s box
            box = obj.getElementBox();
            zSize = diff(box(3,:));
            
            % make rays
            rays = zeros(4, 7);
            rays(:,1:2) = [box(1,1) box(2,1);
                box(1,1) box(2,2);
                box(1,2) box(2,2);
                box(1,2) box(2,1)];
            rays(:,3) = box(3, 2);
            
            % calculate local zMin and zMax to know where is Oz direction
            tmpPoints = GeometryConvertor.getLocalPoint([0 0 box(3,1); 0 0 box(3,2)], obj.CS);
            if tmpPoints(1, 3) > tmpPoints(2, 3)
                zVec = -1;
            else
                zVec = 1;
            end
            
            rays(:,4:6) = repmat(GeometryConvertor.getGlobalVector([0 0 -zVec], obj.CS), [4 1]);
            rays(:,7) = 1;
            
            % check intersections
            [t, ~] = registrator.getInter(rays);
            
            if ~isempty( find(t <=  zSize, 1) )
                res = false;
            else
                res = true;
            end
        end
        
        function surfacesChanged(obj)
            for i = 1 : length(obj.surfaces)
                addlistener(obj.surfaces{i}, 'SurfaceChanged', @obj.surfaceChangedEventListener);
            end
            obj.updateSurfs();
        end
        
        function initListeners(obj)
            if numel(obj.surfacesListeners) > 0
                for i = 1:numel(obj.surfacesListeners)
                    delete(obj.surfacesListeners{i});
                end
                obj.surfacesListeners = {};
            end
            for i = 1:length(obj.surfaces)
                obj.surfacesListeners{i} = addlistener(obj.surfaces{i}, 'SurfaceChanged', @obj.surfaceChangedEventListener);
            end
        end
        function replaceSurface(obj, ind, newSurf)
            obj.surfaces{ind} = newSurf;
            if numel(obj.surfacesListeners) >= ind  && ~isempty(obj.surfacesListeners{ind})
                delete(obj.surfacesListeners{ind});
            end
            obj.surfacesListeners{ind} = addlistener(newSurf, 'SurfaceChanged', @obj.surfaceChangedEventListener);
        end
        
    end % methods
    
    methods (Access = public, Abstract)
        setSurfaces(obj, surfaces);
    end % methods
    
    methods(Static = true, Access = public)
        
        function lensType = getLensTypeByName(name)
            switch name
                case 'LensAxisymTIRAspheric'
                    lensType = 'Lens01';
                case 'LensAxisymTwoAspheric'
                    lensType = 'Lens02';
                case 'LensFreeformNoDome'
                    lensType = 'Lens03';
                case {'LensFreeformSphDome', 'LensFreeformCCO', 'SolidElement'}
                    lensType = 'Lens04';
                case {'LensAxisymTIRUpperFlat'}
                    lensType = 'Lens05';
                case 'LensAxisymTIRFreeform'
                    lensType = 'Lens06';
                case 'LensAxisymFresnelTIR'
                    lensType = 'Lens07';
                case 'MirrorFreeform'
                    lensType = 'Lens08';
                case {'LensTwoFreeform', 'LensAsphericNoDome', 'LensCustom'}
                    lensType = 'Lens09';
                case 'LensNactus'
                    lensType = 'Lens10';
                case {'LensTIRNactus', 'LensFreeformTIRUpperFlat'}
                    lensType = 'Lens11';
                otherwise
                    lensType = 'none';
            end % switch
        end % function
        
    end % methods
    
    methods (Access = protected)
        
        function surfaceChangedEventListener(objSend, eventSrc, eventData)
            try
                objSend.updateSurfs();
            catch ex
                if isa(eventData, 'GeometryEventData')
                    objSend.errorHandle = ex;
                    curSurfs = objSend.surfaces;
                    for i = 1 : length(curSurfs)
                        if (curSurfs{i} == eventSrc)
                            Serializator.update(curSurfs{i}, eventData.bckpStruct);
                        end
                    end
                end
            end
        end
        
        function updateSurfs(~)
        end
        
    end
    
end % classdef

