classdef (Sealed) OpticalScene < ICloneable & ISerializable
    
    properties (Access = {?Serializator})
        solidElements;
        surfaceGroupElements;
        sources;
        registrators;
        envMaterial;
    end % properties
    
    properties (Dependent)
        SolidElements;
        NumSolidElements;
        SurfaceGroupElements;
        NumSurfaceGroupElements;
        OpticalElements;
        NumOpticalElements;
        Sources;
        NumSources;
        Registrators;
        NumRegistrators;
        EnvMaterial;
        CurrentFlux;
    end % properties
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        function value = get.SolidElements(obj)
            value = obj.solidElements;
        end % function

        function value = get.NumSolidElements(obj)
            value = length(obj.solidElements);
        end % function
        
        function value = get.SurfaceGroupElements(obj)
            value = obj.surfaceGroupElements;
        end % function
        
        function value = get.NumSurfaceGroupElements(obj)
            value = length(obj.surfaceGroupElements);
        end % function

        function value = get.OpticalElements(obj)
            if isempty(obj.solidElements)
                if isempty(obj.surfaceGroupElements)
                    value = {};
                else
                    value = obj.surfaceGroupElements;
                end % if
            else
                if isempty(obj.surfaceGroupElements)
                    value = obj.solidElements;
                else
                    value = {obj.solidElements{:}, obj.surfaceGroupElements{:}};
                end % if
            end % if
        end % function
        
        function value = get.NumOpticalElements(obj)
            value = length(obj.solidElements) + length(obj.surfaceGroupElements);
        end % function
        
        function value = get.Sources(obj)
            value = obj.sources;
        end % function
        
        function value = get.NumSources(obj)
            value = length(obj.sources);
        end % function

        function value = get.Registrators(obj)
            value = obj.registrators;
        end % function
        
        function value = get.NumRegistrators(obj)
            value = length(obj.registrators);
        end % function

        function value = get.EnvMaterial(obj)
            value = obj.envMaterial;
        end % function
        
        function value = get.CurrentFlux(obj)
            value = 0;
            for i = 1:obj.NumSources
                if obj.sources{i}.IsEnabled
                    value = value + obj.sources{i}.Flux;
                end
            end
        end
        
        function set.SolidElements(obj, value)
            ValueChecker.checkClass(value, 'SolidElement');
            obj.solidElements = value;
        end % function

        function set.SurfaceGroupElements(obj, value)
            ValueChecker.checkClass(value, 'OpticalElement');
            obj.surfaceGroupElements = value;
        end % function

        function set.OpticalElements(obj, value)
            ValueChecker.checkClass(value, 'OpticalElement');
            % count solids, surfaceGroups
            numSolids = 0; numSurfaceGroups = 0;
            for i = 1 : numel(value)
                if isa(value{i}, 'SolidElement')
                    numSolids = numSolids + 1;
                end % if
                if isa(value{i}, 'SurfaceGroupElement')
                    numSurfaceGroups = numSurfaceGroups + 1;
                end % if
            end % for
            % push solids, surfaceGroups
            solids = cell(numSolids, 1); curSolidNum = 1;
            surfaceGroups = cell(numSurfaceGroups, 1); curSurfaceGroupNum = 1;
            for i = 1 : numel(value)
                if isa(value{i}, 'SolidElement')
                    solids(curSolidNum) = value(i);
                    curSolidNum = curSolidNum + 1;
                end % if
                if isa(value{i}, 'SurfaceGroupElement')
                    surfaceGroups(curSurfaceGroupNum) = value(i);
                    curSurfaceGroupNum = curSurfaceGroupNum + 1;
                end % if
            end % for
            obj.solidElements = solids;
            obj.surfaceGroupElements = surfaceGroups;
        end % function
        
        function set.Sources(obj, value)
            ValueChecker.checkClass(value, 'Source');
            obj.sources = value;
        end % function
        
        function set.Registrators(obj, value)
            ValueChecker.checkClass(value, 'Surface');
            for i = 1:numel(value)
                ValueChecker.checkConstraint(value{i}.IsRegistrator, '==', true);
            end
            obj.registrators = value;
        end % function
                
    end % methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PUBLIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access = public)
        
        function obj = OpticalScene(varArg1)
			if nargin == 1
				if Serializator.isSerializedObj(varArg1)
					ValueChecker.checkDeserializationPossibility(obj, varArg1);
					obj = Serializator.load(obj, varArg1);
				end % if
			else
				obj.envMaterial = Material.Air;
				obj.solidElements = {};
				obj.surfaceGroupElements = {};
			end % if
        end % function
        
        function obj = addOpticalElement(obj, oe)
            ValueChecker.checkClass(oe, 'OpticalElement');
            if isa(oe, 'SolidElement')
                obj.solidElements = [{oe} obj.solidElements];
            end % if
            if isa(oe, 'SurfaceGroupElement')
                obj.surfaceGroupElements{obj.NumSurfaceGroupElements + 1} = oe;
            end % if
        end % function
        
        function obj = addSource(obj, source)
            ValueChecker.checkClass(source, 'Source');
            obj.sources{end+1} = source;
        end % function
        
        function obj = addRegistrator(obj, registrator)
            ValueChecker.checkClass(registrator, 'Surface');
            ValueChecker.checkConstraint(registrator.IsRegistrator, '==', true);
            obj.registrators{end+1} = registrator;
        end % function
        
        function objNew = clone(obj)
            objNew = OpticalScene();
            
            % solids
            for i = 1 : obj.NumSolidElements
                objNew.addOpticalElement(obj.solidElements{i}.clone());
            end % for
            % surf groups
            for i = 1 : obj.NumSurfaceGroupElements
                objNew.addOpticalElement(obj.surfaceGroupElements{i}.clone());
            end % for
            % sources
            for i = 1 : obj.NumSources
                objNew.addSource(obj.sources{i}.clone());
            end % for
            % registrators
            for i = 1 : obj.NumRegistrators
                newReg = obj.registrators{i}.clone();
                newReg.IsRegistrator = obj.registrators{i}.IsRegistrator;
                newReg.Properties = obj.registrators{i}.Properties;
                objNew.addRegistrator(newReg);                
            end % for
            objNew.envMaterial = obj.envMaterial.clone();
        end % function   
        
    end % methods
    
end % classdef

