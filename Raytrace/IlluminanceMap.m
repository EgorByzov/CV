classdef (Sealed) IlluminanceMap < ICloneable & ISerializable
    properties (Access = {?Serializator})
        distr;
        x;
        y;
        symDistr;
        rotSymProf;
        traces;
        registrator;
        options;
    end
    
    properties (Dependent)
        Distr;
        X;
        Y;
        RotationalSymmetryProfile;
        Traces;
        Registrator;
        Options;
    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods 
        
        function value = get.Distr(obj)
            value = obj.symDistr;
        end
        
        function value = get.X(obj)
            value = obj.x;
        end
        
        function value = get.Y(obj)
            value = obj.y;
        end
        
        function value = get.RotationalSymmetryProfile(obj)
            value = obj.rotSymProf;
        end
        
        function value = get.Traces(obj)
            value = obj.traces;
        end % function
        
        function value = get.Registrator(obj)
            value = obj.registrator;
        end % function
        
        function value = get.Options(obj)
            value = obj.options.clone();
        end % function
        
        function set.Options(obj, value)
            ValueChecker.checkClass(value, 'IlluminanceMapOptions');
            ValueChecker.checkAndUpdateSize(value, [1 1]);
            
            if obj.options.NumPointsX ~= value.NumPointsX || ...
                    obj.options.NumPointsY ~= value.NumPointsY || ...
                    obj.options.SigX ~= value.SigX || ...
                    obj.options.SigY ~= value.SigY || ...
                    obj.options.IsSmoothing ~= value.IsSmoothing
                [obj.x, obj.y, obj.distr] = DistributionMapProcessor.makeIlluminanceMap(obj.traces, obj.registrator, value);
                [obj.symDistr, obj.rotSymProf] = IlluminanceMap.applySymmetry(obj.x, obj.y, obj.distr, value.Symmetry);
            elseif obj.options.Symmetry ~= value.Symmetry
                [obj.symDistr, obj.rotSymProf] = IlluminanceMap.applySymmetry(obj.x, obj.y, obj.distr, value.Symmetry);
            end            
            obj.options = value;
        end % function
        
        function set.Traces(obj, value)
            ValueChecker.checkNumerals(value);
            ValueChecker.checkAndUpdateSize(value, [Inf 3]);
            
            if ~isequal(value, obj.traces)
                obj.traces = value;
                [obj.x, obj.y, obj.distr] = DistributionMapProcessor.makeIlluminanceMap(value, obj.registrator, obj.options);
                [obj.symDistr, obj.rotSymProf] = IlluminanceMap.applySymmetry(obj.x, obj.y, obj.distr, obj.options.Symmetry);
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PUBLIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access = public)        
        % constructor
        function obj = IlluminanceMap(varArg1, registrator, options, x, y, distr)
			% varArg1 - traces
			if nargin == 1 && Serializator.isSerializedObj(varArg1)
				ValueChecker.checkDeserializationPossibility(obj, varArg1);
				obj = Serializator.load(obj, varArg1);
            else
                if ~isempty(varArg1)
                    ValueChecker.checkNumerals(varArg1);
                    ValueChecker.checkAndUpdateSize(varArg1, [Inf 3]);
                end
				ValueChecker.checkClass(registrator, 'Surface');
				ValueChecker.checkClass(options, 'IlluminanceMapOptions');
				
				ValueChecker.checkAndUpdateSize(registrator, [1 1]);
				ValueChecker.checkAndUpdateSize(options, [1 1]);

				ValueChecker.checkConstraint(registrator.IsRegistrator, '==', true);

				obj.traces = varArg1;
				obj.registrator = registrator;
				obj.options = options;

				if nargin == 3
					[obj.x, obj.y, obj.distr] = DistributionMapProcessor.makeIlluminanceMap(varArg1, registrator, options);
				elseif nargin == 6
					ValueChecker.checkNumerals(x);
					ValueChecker.checkNumerals(y);
					ValueChecker.checkNumerals(distr);
					ValueChecker.checkConstraint(distr, '>=', 0);

					ValueChecker.checkConstraint(size(distr), '==', size(x));
					ValueChecker.checkConstraint(size(distr), '==', size(y));

					obj.x = x;
					obj.y = y;
					obj.distr = distr;
				end % if
				[obj.symDistr, obj.rotSymProf] = IlluminanceMap.applySymmetry(obj.x, obj.y, obj.distr, options.Symmetry);
			end % if
        end % function
          
        function objNew = clone(obj)
            objNew = IlluminanceMap(obj.traces, obj.registrator, obj.options.clone(), obj.x, obj.y, obj.distr);
        end % function

    end % methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PRIVATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Static = true, Access = private)
        function [symDistr, rotSymFunc] = applySymmetry(x, y, distr, symmetry)            
            switch (symmetry)
                case Constants.IlluminanceMapSymmetries.QUADRANT
                    symDistr = distr + flipud(distr);
                    symDistr = (symDistr + fliplr(symDistr))/4;
                    rotSymFunc = [];
                case Constants.IlluminanceMapSymmetries.UD
                    symDistr = (distr + flipud(distr))/2;
                    rotSymFunc = [];
                case Constants.IlluminanceMapSymmetries.LR
                    symDistr = (distr + fliplr(distr))/2;
                    rotSymFunc = [];
                case Constants.IlluminanceMapSymmetries.ROTATIONAL
                    % prepare params
                    irradR = sqrt(x.^2 + y.^2);
                    numPoints = ceil(sqrt(size(x,1)^2 + size(x,2)^2)/2);
                    Rmax = max(max(irradR));
                    dR = Rmax / (numPoints - 1);
                    Rnodes = 0 : dR : Rmax;
                    newIllum = zeros(size(Rnodes));
                    goodInds = true(size(newIllum));
                    
                    % fill irradSpline
                    for i = 1:numPoints
                        inds = find( (irradR > Rnodes(i) - dR/2) & (irradR < Rnodes(i) + dR/2) );
                        if ~isempty(inds)
                            newIllum(i) = sum(sum(distr(inds))) / numel(inds);
                        else
                            goodInds(i) = false;
                        end
                    end
                    newIllum = newIllum(goodInds);
                    Rnodes = Rnodes(goodInds);
                    if Rnodes(1) ~= 0
                        Rnodes = [0 Rnodes];
                        newIllum = [newIllum(1) newIllum];
                    end
                    if Rnodes(end) < Rmax
                        Rnodes = [Rnodes Rmax];
                        newIllum = [newIllum newIllum(end)];
                    end
                    rotSymFunc = @(x, y) interp1(Rnodes, newIllum, sqrt(x.^2 + y.^2));
                    symDistr = rotSymFunc(x, y);
                otherwise
                    symDistr = distr;
                    rotSymFunc = [];
            end
        end
    end
end