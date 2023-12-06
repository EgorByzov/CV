classdef SolidElement < OpticalElement & IDirectRaytraceable
    
    properties (Access = {?SolidElement, ?Serializator})      
        innMaterial;
    end % properties
        
    properties (Dependent)
        InnMaterial;
    end % properties
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods 
        
        function value = get.InnMaterial(obj)
            value = obj.innMaterial;
        end % function
        
        function set.InnMaterial(obj, innMaterial)
            ValueChecker.checkClass(innMaterial, 'Material');
            obj.innMaterial = innMaterial;
        end % function
        
    end % methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PUBLIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    methods (Access = public)
    
        function obj = SolidElement(varArg1, innMat)
			% varArg1 - surface
			if nargin == 1
				if Serializator.isSerializedObj(varArg1)
					ValueChecker.checkDeserializationPossibility(obj, varArg1);
					obj = Serializator.load(obj, varArg1);
				end
			elseif nargin == 2
                ValueChecker.checkClass(varArg1, 'Surface');
                ValueChecker.checkClass(innMat, 'Material');
                obj.innMaterial = innMat;
                obj.surfaces = varArg1;
                obj.cs = CoordinateSystem.Global;
			end % if
        end % function
                    
        function objNew = clone(obj)
            newSurfaces = obj.surfaces;
            numSurfs = numel(newSurfaces);
            for i = 1 : numSurfs
                newSurfaces{i} = newSurfaces{i}.clone();
            end % for
            objNew = SolidElement(newSurfaces, obj.innMaterial.clone());
            objNew.cs = obj.cs.clone();
        end % function        
    
        function setSurfaces(obj, surfaces)
            obj.surfaces = surfaces;
            obj.surfacesChanged();
        end
        
    end % methods
	
end % classdef