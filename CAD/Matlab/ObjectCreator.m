classdef (Sealed) ObjectCreator < handle & ISerializable
    
    properties (Access = {?Serializator})
        % Optical Scene objects
        sources = {}; % {source params}
        lenses = {}; % {className ZAstruct computedResult}
        registrators = {};
        
        % Optical Task Objects
        requiredDistributions = {}; % {reqDistr params}
        
        % inds
        curSource = [];
        curLense = [];
        curRegistrator = [];
        curReqDistr = [];        
    end
    
    properties (Dependent)
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PUBLIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    methods (Access = public)
        
		function obj = ObjectCreator(varArg1)
			if nargin == 1
				if Serializator.isSerializedObj(varArg1)
					ValueChecker.checkDeserializationPossibility(obj, varArg1);
					obj = Serializator.load(obj, varArg1);
				end % if
			end
		end
		
        function [source, lense, requiredDistribution, rays, map] = startStandardCreationWizard(obj, parentHandle, project, properties)
            % initialize vars
            step = 1;
            set(0, 'units', 'pixels');
            screensize = get(0, 'screensize');
            x = (screensize(3) - 560)/2;
            y = (screensize(4) - 460)/2;
            setappdata(parentHandle, 'childPosition', [x y 560 460]);
            % init default values
            for i = 1 : numel(obj.registrators)
                obj.registrators{i}.CS.Center = obj.registrators{i}.CS.Center .* [0 0 1];
                obj.registrators{i}.CS.setAxes([0 1 0], [0 0 1]); 
            end

            % start dialog
            status = obj.runCurrentStandardWizardStep(step, project, parentHandle, properties);            
            % main loop
            while ~strcmp(status, 'Finish') && ~strcmp(status, 'Cancel')                
                switch status                    
                    case 'Back'                        
                        step = step - 1;
                        status = obj.runCurrentStandardWizardStep(step, project, parentHandle, properties);                        
                    case 'Next'                        
                        step = step + 1;
                        status = obj.runCurrentStandardWizardStep(step, project, parentHandle, properties);                        
                end % switch                
            end % while
            
            % delete trash
            if isappdata(parentHandle, 'plots')
                deleteHandles(getappdata(parentHandle, 'plots'));
                rmappdata(parentHandle, 'plots');
            end 
            rmappdata(parentHandle, 'childPosition');
            
            if strcmp(status, 'Finish')
                results = obj.lenses{obj.curLense, 3};
                source = results.source;
                lense = results.lens;
                requiredDistribution = results.reqDistr;
                rays = results.rays;
                map = results.map;
                if isa(requiredDistribution, 'RequiredLightDistribution')
                    if isa(requiredDistribution, 'RequiredIlluminanceDistribution')
                        perfRaytr = requiredDistribution.makeRaytracePerformance(source.Flux, map.Traces);
                    else
                        perfRaytr = requiredDistribution.makeRaytracePerformance(source.Flux, results.rays);
                    end
                    if perfRaytr.effShaping > 0
                        requiredDistribution.RequiredEff = perfRaytr.effShaping/100;
                    end
                end
            else
                source = [];
                lense = [];
                requiredDistribution = [];
                rays = [];
                map = [];
            end
        end
       
        % % % % % % % % % % EDIT % % % % % % % % % % % % % % % % % % % % %%
        
        function reqDistr = editRequiredDistribution(obj, reqDistr, parentHandle, project, properties)
            [~, num] = obj.addReqDistr(reqDistr);
            params = obj.requiredDistributions{num, 2};
            flux = project.OS.CurrentFlux;
            if project.AreRaysTraced
                rays = project.TracedRays;   
                if isa(reqDistr, 'RequiredIlluminanceDistribution')
                    pSet = reqDistr.Registrator.getPrimitiveSet(project.ApproximOptions);
                    traces = Raytracer.makeTraces(pSet, rays);
                else
                    traces = [];
                end
            else
                rays = [];
                traces = [];
            end
            
            if isa(reqDistr, 'RequiredLightDistribution')
                switch class(reqDistr.LightDistribution)
                    case {'IlluminanceAxisym', 'Illuminance1ProfileInComplexArea'}
                        frmFunc = @frmEditStandartReqIrradModel;
                        raytraceResult = traces;
                    case 'IntensityAxisym'
                        frmFunc = @frmEditStandardReqIntenModel;
                        raytraceResult = rays;
                    case 'Illuminance2ProfilesInSimpleArea'
                        frmFunc = @frmEdit2ProfileReqIllum;
                        raytraceResult = traces;
                    case 'Intensity2Profiles'
                        frmFunc = @frmEdit2ProfileReqInten;
                        raytraceResult = rays;
                    case 'IlluminanceInterp3D'
                        frmFunc = @frmEditInterp3DReqIrradModel;
                        raytraceResult = traces;
                    case 'IntensityInterpInPlaneDomain'
                        frmFunc = @frmEditIESOnRoadReqIrradModel;
                        raytraceResult = rays;
                    otherwise
                        frmFunc = @frmEditReqDistrModel;
                        raytraceResult = rays;
                end
            elseif isa(reqDistr, 'RequiredCollimatedDistribution')
                frmFunc = @frmEditCollimatedDistrModel;
                raytraceResult = rays;
            end 
            ch = frmFunc('frmParent', parentHandle,...
                         'reqDistr', reqDistr,...
                         'params', params,...
                         'mode', 'Edit',...
                         'properties', properties,...
                         'raytraceResult', raytraceResult,...
                         'flux', flux);
            uiwait(ch);
            status = getappdata(parentHandle, 'childResult');
            objRes = getappdata(parentHandle, 'childObject');
            obj.requiredDistributions{num, 2} = objRes;     
            if ~reqDistr.equals(obj.requiredDistributions{num, 1})
                obj.addReqDistr(reqDistr);
            end
            
            if strcmp(status, 'Redefine')
                oldReqDistr =  reqDistr;
                while ~strcmp(status, 'Cancel') && ~strcmp(status, 'Finish') && ~strcmp(status, 'Ok')
                    switch status
                        case {'Redefine', 'Back'}
                            if isa(reqDistr, 'RequiredLightDistribution')
                                reqDistrType = ['LightDistribution:' class(reqDistr.LightDistribution)];
                            elseif isa(reqDistr, 'RequiredCollimatedDistribution')
                                reqDistrType = 'CollimatingDistribution';
                            end
                            ch = frmChooseReqIrradModelType('frmParent', parentHandle,...
                                'reqDistrType', reqDistrType, 'mode', 'Redefine', ...
                                'properties', properties);
                            uiwait(ch);
                            objRes = getappdata(parentHandle, 'childObject');
                        case {'Next'}
                            [reqDistr, num] = obj.addReqDistr(objRes);
                            params = obj.requiredDistributions{num, 2};
                            if isa(reqDistr, 'RequiredLightDistribution')
                                switch class(reqDistr.LightDistribution)
                                    case {'IlluminanceAxisym', 'Illuminance1ProfileInComplexArea'}
                                        frmFunc = @frmEditStandartReqIrradModel;
                                        
                                        if ~isempty(rays) && isempty(traces)
                                            pSet = reqDistr.Registrator.getPrimitiveSet(project.ApproximOptions);
                                            traces = Raytracer.makeTraces(pSet, rays);
                                        end
                                        raytraceResult = traces;
                                    case 'IntensityAxisym'
                                        frmFunc = @frmEditStandardReqIntenModel;
                                        raytraceResult = rays;
                                    case 'Illuminance2ProfilesInSimpleArea'
                                        frmFunc = @frmEdit2ProfileReqIllum;
                                        if ~isempty(rays) && isempty(traces)
                                            pSet = reqDistr.Registrator.getPrimitiveSet(project.ApproximOptions);
                                            traces = Raytracer.makeTraces(pSet, rays);
                                        end
                                        raytraceResult = traces;
                                    case 'Intensity2Profiles'
                                        frmFunc = @frmEdit2ProfileReqInten;
                                        raytraceResult = rays;
                                    case 'IlluminanceInterp3D'
                                        frmFunc = @frmEditInterp3DReqIrradModel;
                                        if ~isempty(rays) && isempty(traces)
                                            pSet = reqDistr.Registrator.getPrimitiveSet(project.ApproximOptions);
                                            traces = Raytracer.makeTraces(pSet, rays);
                                        end
                                        raytraceResult = traces;
                                    case 'IntensityInterpInPlaneDomain'
                                        frmFunc = @frmEditIESOnRoadReqIrradModel;
                                        raytraceResult = rays;
                                end
                            elseif isa(reqDistr, 'RequiredCollimatedDistribution')
                                frmFunc = @frmEditCollimatedDistrModel;
                                raytraceResult = rays;
                            end
                            ch = frmFunc('frmParent', parentHandle,...
                                        'reqDistr', reqDistr,...
                                        'params', params,...
                                        'mode', 'Edit',...
                                        'properties', properties,...
                                        'raytraceResult', raytraceResult,...
                                        'flux', flux);
                            uiwait(ch);
                            objRes = getappdata(parentHandle, 'childObject');
                            obj.requiredDistributions{num, 2} = objRes;
                    end
                    status = getappdata(parentHandle, 'childResult');
                end
                if strcmp(status, 'Cancel')
                    reqDistr = oldReqDistr;
                else
                    reqDistr = reqDistr.clone();
                    if isa(reqDistr, 'RequiredIlluminanceDistribution') && isa(oldReqDistr, 'RequiredIlluminanceDistribution')
                        reqDistr.Registrator = oldReqDistr.Registrator;
                    end
                end
            end
        end
        
        function editRegistrator(obj, registrator, parentHandle, properties)
            [~, num] = obj.addRegistrator(registrator);
            
            switch class(registrator)
                case 'PlanarQuadrangle'
                    frmFunc = @frmEditExitPlane;                   
            end
            ch = frmFunc('frmParent', parentHandle,...
                         'surface', registrator,...
                         'properties', properties);
            uiwait(ch);
            if ~isequal(registrator, obj.registrators{num})
                obj.addRegistrator(registrator);
            end
        end
        
        function source = editSource(obj, source, solid, parentHandle, properties)
            [~, num] = obj.addSource(source);
            params = obj.sources{num, 2};
            
            switch class(source)
                case 'StandardSource'
                    frmFunc = @frmEditStandartLEDModel;
                case 'RayFileSource'
                    frmFunc = @frmEditTPRayFileLEDModel;                    
            end
            ch = frmFunc('frmParent', parentHandle,...
                         'source', source,...
                         'params', params,...
                         'mode', 'Edit',...
                         'lens', solid,...
                         'properties', properties);
            uiwait(ch);
            status = getappdata(parentHandle, 'childResult');
            objRes = getappdata(parentHandle, 'childObject');
            obj.sources{num, 2} = objRes;     
            if ~isequal(source, obj.sources{num, 1})
                obj.addSource(source);
            end
            
            if strcmp(status, 'Redefine source')
                oldSource =  source;
                while ~strcmp(status, 'Cancel') && ~strcmp(status, 'Finish') && ~strcmp(status, 'Ok')
                    switch status
                        case {'Redefine source', 'Back'}
                            ch = frmChooseLEDModelType('frmParent', parentHandle,...
                                                       'sourceType', class(source),...
                                                       'mode', 'Redefine', ...
                                                       'properties', properties);
                            uiwait(ch);
                            objRes = getappdata(parentHandle, 'childObject');
                        case {'Next'}
                            [source, num] = obj.addSource(objRes);
                            params = obj.sources{num, 2};
                            switch objRes
                                case 'StandardSource'
                                    frmFunc = @frmEditStandartLEDModel;
                                case 'RayFileSource'
                                    frmFunc = @frmEditTPRayFileLEDModel;
                            end
                            ch = frmFunc('frmParent', parentHandle,...
                                         'source', source,...
                                         'params', params,...
                                         'mode', 'Redefine source',...
                                         'lens', solid,...
                                         'properties', properties);
                            uiwait(ch);
                            objRes = getappdata(parentHandle, 'childObject');
                            obj.sources{num, 2} = objRes;
                    end
                    status = getappdata(parentHandle, 'childResult');
                end
                if strcmp(status, 'Cancel')
                    source = oldSource;
                else
                    source = source.clone();
                end
            end
        end
        
        function editSolidElement(~, os, parentHandle, properties)
            solid = os.OpticalElements{1};
            if solid.isAxisymmetrical()
                ch = frmSolidEditor2D('frmParent', parentHandle,...
                                      'opticalScene', os,...
                                      'properties', properties);
            else
                ch = frmOEEditor3D('frmParent', parentHandle,...
                                   'opticalScene', os,...
                                   'properties', properties);
            end
            uiwait(ch);
        end
        
        function [solid, rays, map]  = computeOpticalElement(obj, varArg, reqDistr, source, parentHandle, properties)
            raytrOpt = DirectRaytracingOptions.Default;
            if isa(reqDistr, 'RequiredIlluminanceDistribution')
                mapOpt = IlluminanceMapOptions.Default;
            else
                mapOpt = IntensityMapOptions.Default;
            end
            
            % check cur lens
            [oeType, pos] = obj.addOpticalElement(varArg);
            
            if ~isempty(obj.lenses{pos, 3})
                params = obj.lenses{pos, 3};
                if isequal(params.reqDistr, reqDistr) && ...
                        isequal(params.source, source)
                    solidParam = params.lens;
                else
                    solidParam = obj.lenses{pos, 2};
                end
            else
                solidParam = obj.lenses{pos, 2};
            end
            
            s = sprintf('frm%sProperties', oeType);
            func = str2func(s);
            ch = func(...
                        'Position', getappdata(parentHandle, 'childPosition'),...
                        'frmParent', parentHandle,...
                        'solid', solidParam,...
                        'reqDistr', reqDistr,...
                        'source', source,...
                        'raytrOpt', raytrOpt,...
                        'mapOptions', mapOpt,...
                        'mode', 'Recompute',...
                        'properties', properties);
            uiwait(ch);
            
            % get results and process them
            status = getappdata(parentHandle, 'childResult');
            objRes = getappdata(parentHandle, 'childObject');
            
            if isfield(objRes, 'lens')
                if ~strcmp(status, 'Cancel')
                    solid = objRes.lens;
                    rays = objRes.rays;
                    map = objRes.map;
                else
                    solid = [];
                    rays = [];
                    map = [];
                end
                objRes.reqDistr = reqDistr.clone();
                if isa(objRes.reqDistr, 'RequiredIlluminanceDistribution') && ...
                                ~isempty(objRes.reqDistr.Registrator)
                    objRes.reqDistr.Registrator = objRes.reqDistr.Registrator.clone();
                end
                objRes.source = source.clone();
                obj.lenses{obj.curLense, 2} = objRes.lens.ZA;
                obj.lenses{obj.curLense, 3} = objRes;
            else
                solid = [];
                rays = [];
                map = [];
                if ~isequal(obj.lenses{obj.curLense, 2}, objRes)
                    obj.lenses{obj.curLense, 2} = objRes;
                    obj.lenses{obj.curLense, 3} = [];
                end
            end
        end
        
        % % % % % % % % % % ADD % % % % % % % % % % % % % % % % % % % % % %
        
        function [source, pos] = addSource(obj, varArg)
            if isa(varArg, 'Source')
                source = varArg.clone();
                [isSource, pos] = obj.isSuchSourceIn(class(source));
                if isSource            
                    obj.sources{pos, 1} = source;
                    obj.curSource = pos;
                else
                    switch class(varArg)
                        case 'StandardSource'
                            params.spline = [];
                            params.sigma  = 0.4;
                        case 'RayFileSource'
                            params = [];
                    end % switch
                    obj.sources = [obj.sources ; {source params}];
                    obj.curSource = size(obj.sources, 1);
                end
            else
                [isSource, pos] = obj.isSuchSourceIn(varArg);
                if isSource
                    obj.curSource = pos;
                    source = obj.sources{pos, 1};
                else
                    switch varArg
                        case 'StandardSource'
                            source = StandardSource.Default;
                            params.spline = [];
                            params.sigma  = 0.4;
                        case 'RayFileSource'
                            source = RayFileSource.Default;
                            params = [];
                    end % switch
                    obj.sources = [obj.sources ; {source params}];
                    obj.curSource = size(obj.sources, 1);
                end % if
            end % if
            pos = obj.curSource;
        end % function
                
        function [registrator, pos] = addRegistrator(obj, varArg)
            if isa(varArg, 'Surface')
                registrator = varArg.clone();
                                
                [isRegistrator, pos] = obj.isSuchRegistratorIn(class(registrator));
                if isRegistrator
                    for i = 1:size(obj.requiredDistributions, 1)
                        if isa(obj.requiredDistributions{i}, 'RequiredIlluminanceDistribution') && obj.requiredDistributions{i}.Registrator == obj.registrators{pos}
                            obj.requiredDistributions{i}.Registrator = registrator;
                        end
                    end
                    obj.registrators{pos} = registrator;
                    obj.curRegistrator = pos;
                else
                    obj.registrators = [obj.registrators; {registrator}];
                    obj.curRegistrator = size(obj.registrators, 1);
                end
            else
                [isRegistrator, pos] = obj.isSuchRegistratorIn(varArg);
                if isRegistrator
                    registrator = obj.registrators{pos};
                    obj.curRegistrator = pos;
                else
                    switch varArg
                        case 'PlanarQuadrangle'
                            registrator = PlanarQuadrangle.Default;
                            registrator.IsRegistrator = true;
                            obj.registrators = [obj.registrators; {registrator}];
                            obj.curRegistrator = size(obj.registrators, 1);
                    end
                end                
            end
            pos = obj.curRegistrator;
        end
                
        function [reqDistr, pos] = addReqDistr(obj, varArg)
            [isReqDistr, pos] = obj.isSuchReqDistrIn(varArg);
            if isa(varArg, 'RequiredDistribution')
                reqDistr = varArg.clone();
                if isa(reqDistr, 'RequiredIlluminanceDistribution')
                    reqDistr.Registrator = obj.addRegistrator(reqDistr.Registrator);                                  
                end                
                if isReqDistr            
                    obj.requiredDistributions{pos, 1} = reqDistr;
                    obj.curReqDistr = pos;
                else
                    if isa(reqDistr, 'RequiredLightDistribution')
                        switch class(reqDistr.LightDistribution)
                            case 'IlluminanceAxisym'
                                params.spline = [];
                                params.sigma = 100;
                            case 'Illuminance1ProfileInComplexArea'
                                params.spline = [];
                                params.sigma = 100;
                            case 'Illuminance2ProfilesInSimpleArea'
                                params.splineX = [];
                                params.sigmaX = 100;
                                params.splineY = [];
                                params.sigmaY = 100;
                            case 'IntensityAxisym'
                                params.spline = [];
                                params.sigma = 0.4;
                            case 'Intensity2Profiles'
                                params.splineX = [];
                                params.sigmaX = 0.4;
                                params.splineY = [];
                                params.sigmaY = 0.4;
                            case 'IntensityInterpInPlaneDomain'
                                params.spline = Constant(1);                                
                            otherwise
                                params = [];
                        end
                    elseif isa(reqDistr, 'RequiredCollimatedDistribution')
                        params = [];
                    else
                        params = [];
                    end
                    obj.requiredDistributions = [obj.requiredDistributions ; {reqDistr params}];
                    obj.curReqDistr = size(obj.requiredDistributions, 1);
                end
            else
                if isReqDistr
                    obj.curReqDistr = pos;
                    reqDistr = obj.requiredDistributions{obj.curReqDistr, 1};
                    if  isa(reqDistr, 'RequiredIlluminanceDistribution')
                        [isRegistrator, pos] = obj.isSuchRegistratorIn(class(reqDistr.Registrator));
                        if isRegistrator
                            obj.curRegistrator = pos;
                        else
                            obj.registrators = [obj.registrators; {reqDistr.Registrator}];
                            obj.curRegistrator = size(obj.registrators, 1);
                        end
                    end
                else
                    switch varArg
                        case 'LightDistribution:IlluminanceAxisym'
                            params.spline = [];
                            params.sigma = 100;
                            newLightDistr = IlluminanceAxisym.Default;
                            newRegistrator = obj.addRegistrator('PlanarQuadrangle');
                            reqDistr = RequiredIlluminanceDistribution(newLightDistr, newRegistrator);                        
                        case 'LightDistribution:Illuminance1ProfileInComplexArea'
                            params.spline = [];
                            params.sigma = 100;
                            newLightDistr = Illuminance1ProfileInComplexArea.Default;
                            newRegistrator = obj.addRegistrator('PlanarQuadrangle');
                            reqDistr = RequiredIlluminanceDistribution(newLightDistr, newRegistrator);
                        case 'LightDistribution:Illuminance2ProfilesInSimpleArea'
                            params.splineX = [];
                            params.sigmaX = 100;
                            params.splineY = [];
                            params.sigmaY = 100;
                            newLightDistr = Illuminance2ProfilesInSimpleArea.Default;
                            newRegistrator = obj.addRegistrator('PlanarQuadrangle');
                            reqDistr = RequiredIlluminanceDistribution(newLightDistr, newRegistrator);
                        case 'LightDistribution:IntensityAxisym'
                            params.spline = [];
                            params.sigma = 0.4;
                            newLightDistr = IntensityAxisym.Default;
                            reqDistr = RequiredIntensityDistribution(newLightDistr);
                        case 'LightDistribution:Intensity2Profiles'
                            params.splineX = [];
                            params.sigmaX = 0.4;
                            params.splineY = [];
                            params.sigmaY = 0.4;
                            newLightDistr = Intensity2Profiles.Default;
                            reqDistr = RequiredIntensityDistribution(newLightDistr);
                        case 'LightDistribution:IlluminanceInterp3D'
                            params = [];
                            newLightDistr = IlluminanceInterp3D.Default;
                            newRegistrator = obj.addRegistrator('PlanarQuadrangle');
                            reqDistr = RequiredIlluminanceDistribution(newLightDistr, newRegistrator);
                        case 'LightDistribution:IntensityByIESOnRoad'
                            params.spline = Constant(1);
                            newLightDistr = IntensityInterpInPlaneDomain.Default;
                            reqDistr = RequiredIntensityByIllumDistribution(newLightDistr, 1);
                        case 'CollimatingDistribution'
                            reqDistr = RequiredCollimatedDistribution();
                            params = [];
                    end
                    
                    obj.requiredDistributions = [obj.requiredDistributions ; {reqDistr params}];
                    obj.curReqDistr = size(obj.requiredDistributions, 1);
                end
            end
            pos = obj.curReqDistr;
        end
       
        function [opticalElement, pos] = addOpticalElement(obj, varArg)
            if isa(varArg, 'OpticalElement')
                opticalElement = class(varArg);
                [isLens, pos] = obj.isSuchLenseIn(opticalElement);
                if isLens            
                    obj.lenses{pos, 2} = varArg.ZA;
                    obj.lenses{pos, 3} = [];
                    obj.curLense = pos;
                else
                    obj.lenses = [obj.lenses; {opticalElement varArg.ZA []}];
                    obj.curLense = size(obj.lenses, 1);
                end
            else                
                opticalElement = varArg;
                [isLens, pos] = obj.isSuchLenseIn(varArg);
                if isLens
                    obj.curLense = pos;
                else
                    lensFunc = str2func([opticalElement '.DefaultZA']);
                    obj.lenses = [obj.lenses; {opticalElement lensFunc() []}];
                    obj.curLense = size(obj.lenses, 1);
                end % if
            end % if
            pos = obj.curLense;
        end % function
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PRIVATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access = private)
        
        function [isSource, pos] = isSuchSourceIn(obj, className)
            isSource = false;
            pos = 0;
            
            for i = 1:size(obj.sources, 1)
                if isa(obj.sources{i, 1}, className)
                    isSource = true;
                    pos = i;
                    return;
                end
            end            
        end
        
        function [isLense, pos] = isSuchLenseIn(obj, className)
            isLense = false;
            pos = 0;
            
            for i = 1:size(obj.lenses, 1)
                if strcmp(obj.lenses{i, 1}, className)
                    isLense = true;
                    pos = i;
                    return;
                end
            end
        end
        
        function [isRegistrator, pos] = isSuchRegistratorIn(obj, className)
            isRegistrator = false;
            pos = 0;
            
            for i = 1:size(obj.registrators, 1)
                if isa(obj.registrators{i}, className)
                    isRegistrator = true;
                    pos = i;
                    return;
                end
            end
        end
                        
        function [isReqDistr, pos] = isSuchReqDistrIn(obj, className)
            if ~ischar(className)
                if isa(className, 'RequiredCollimatedDistribution')
                    className = 'CollimatingDistribution';
                elseif isa(className, 'RequiredLightDistribution')
                    className = ['LightDistribution:' class(className.LightDistribution)];
                end
            end
            isReqDistr = false;
            pos = 0;
            
            for i = 1:size(obj.requiredDistributions, 1)
                if isa(obj.requiredDistributions{i, 1}, 'RequiredLightDistribution')
                    if strcmp(className, ['LightDistribution:' class(obj.requiredDistributions{i, 1}.LightDistribution)])
                        isReqDistr = true;
                        pos = i;
                        return;
                    end
                elseif isa(obj.requiredDistributions{i, 1}, 'RequiredCollimatedDistribution')...
                            && strcmp(className, 'CollimatingDistribution')
                    isReqDistr = true;
                    pos = i;
                    return;
                end
            end
        end
                
        function status = runCurrentStandardWizardStep(obj, step, project, parentHandle, properties)
            
            switch step
                case 1
                    if ~isempty(obj.curSource)
                        sourceType = class(obj.sources{obj.curSource});
                    else
                        sourceType = [];
                    end
                    ch = frmChooseLEDModelType(...
                        'Position', getappdata(parentHandle, 'childPosition'),...
                        'frmParent', parentHandle,...
                        'sourceType', sourceType, 'mode', 'Project creation',...
                        'properties', properties);
                case 2
                    source = obj.sources{obj.curSource, 1};
                    source.EmitSol.CS = CoordinateSystem.Global;
                    params = obj.sources{obj.curSource, 2};
                    switch class(source)
                        case 'StandardSource'
                            ch = frmEditStandartLEDModel(...
                                'Position', getappdata(parentHandle, 'childPosition'), ...
                                'source', source, 'params', params,...
                                'mode', 'Project creation',...
                                'frmParent', parentHandle,...
                                'properties', properties);
                        case 'RayFileSource'
                            ch = frmEditTPRayFileLEDModel(...
                                'Position', getappdata(parentHandle, 'childPosition'),...
                                'frmParent', parentHandle,...
                                'source', source,...
                                'mode', 'Project creation',...
                                'properties', properties);
                    end
                case 3
                    if ~isempty(obj.curReqDistr)
                        reqDistr = obj.requiredDistributions{obj.curReqDistr, 1};
                        if isa(reqDistr, 'RequiredLightDistribution')
                            reqDistrType = ['LightDistribution:' class(reqDistr.LightDistribution)];
                        elseif isa(reqDistr, 'RequiredCollimatedDistribution')
                            reqDistrType = 'CollimatingDistribution';
                        end
                    else
                        reqDistrType = [];
                    end
                    ch = frmChooseReqIrradModelType(...
                        'Position', getappdata(parentHandle, 'childPosition'),...
                        'frmParent', parentHandle,...
                        'reqDistrType', reqDistrType, 'mode', 'Project creation',...
                        'properties', properties);
                case 4
                    reqDistr = obj.requiredDistributions{obj.curReqDistr, 1};  
                    params = obj.requiredDistributions{obj.curReqDistr, 2};
                    
                    if isa(reqDistr, 'RequiredLightDistribution')
                        reqDistr.RequiredEff = 1;
                        switch class(reqDistr.LightDistribution)
                            case {'IlluminanceAxisym', 'Illuminance1ProfileInComplexArea'}
                                ch = frmEditStandartReqIrradModel(...
                                    'Position', getappdata(parentHandle, 'childPosition'),...
                                    'frmParent', parentHandle,...
                                    'reqDistr', reqDistr,...
                                    'params', params,...
                                    'mode', 'Project creation',...
                                    'properties', properties);
                            case 'IntensityAxisym'
                                ch = frmEditStandardReqIntenModel(...
                                    'Position', getappdata(parentHandle, 'childPosition'),...
                                    'frmParent', parentHandle,...
                                    'reqDistr', reqDistr,...
                                    'params', params,...
                                    'mode', 'Project creation',...
                                    'properties', properties);
                            case 'Illuminance2ProfilesInSimpleArea'
                                ch = frmEdit2ProfileReqIllum(...
                                    'Position', getappdata(parentHandle, 'childPosition'),...
                                    'frmParent', parentHandle,...
                                    'reqDistr', reqDistr,...
                                    'params', params,...
                                    'mode', 'Project creation',...
                                    'properties', properties);
                            case 'Intensity2Profiles'
                                ch = frmEdit2ProfileReqInten(...
                                    'Position', getappdata(parentHandle, 'childPosition'),...
                                    'frmParent', parentHandle,...
                                    'reqDistr', reqDistr,...
                                    'params', params,...
                                    'mode', 'Project creation',...
                                    'properties', properties);
                            case 'IlluminanceInterp3D'
                                ch = frmEditInterp3DReqIrradModel(...
                                    'Position', getappdata(parentHandle, 'childPosition'),...
                                    'frmParent', parentHandle,...
                                    'reqDistr', reqDistr,...
                                    'params', params,...
                                    'mode', 'Project creation',...
                                    'properties', properties);
                            case 'IntensityInterpInPlaneDomain'
                                ch = frmEditIESOnRoadReqIrradModel(...
                                    'Position', getappdata(parentHandle, 'childPosition'),...
                                    'frmParent', parentHandle,...
                                    'reqDistr', reqDistr,...
                                    'params', params,...
                                    'mode', 'Project creation',...
                                    'properties', properties);
                        end
                    elseif isa(reqDistr, 'RequiredCollimatedDistribution')
                         ch = frmEditCollimatedDistrModel(...
                            'Position', getappdata(parentHandle, 'childPosition'),...
                            'frmParent', parentHandle,...
                            'reqDistr', reqDistr,...
                            'params', params,...
                            'mode', 'Project creation',...
                            'properties', properties);                       
                    end
                case 5
                    if ~isempty(obj.curLense)
                        lensType = obj.lenses{obj.curLense, 1};
                    else
                        lensType = [];
                    end  
                    reqDistr = obj.requiredDistributions{obj.curReqDistr, 1}.clone();
                    ch = frmChooseLensType(...
                        'Position', getappdata(parentHandle, 'childPosition'),...
                        'frmStart', parentHandle,...
                        'lensType', lensType,...
                        'reqDistr', reqDistr,...
                        'properties', properties);
                case 6
                    reqDistr = obj.requiredDistributions{obj.curReqDistr, 1}.clone();                    
                    source = obj.sources{obj.curSource, 1}.clone;    
                    raytrOpt = project.DirectRaytrOptions;
                    
                    s = sprintf('frm%sProperties', obj.lenses{obj.curLense, 1});
                    func = str2func(s);
                    
                    if isa(reqDistr, 'RequiredIlluminanceDistribution')
                        registrator = reqDistr.Registrator;
                        mapOptions = project.IllumMapOptions;
                        lightDistr = reqDistr.LightDistribution;
                        switch class(lightDistr)
                            case {'IlluminanceAxisym', 'Illuminance1ProfileInComplexArea'}
                                RMax = lightDistr.Profile.RMax;
                                if isa(registrator, 'PlanarQuadrangle')
                                    registrator.Width = RMax * 2.8;
                                    registrator.Height = RMax * 2.8;
                                end
                            case 'Illuminance2ProfilesInSimpleArea'
                                RMax = max(lightDistr.ProfileX.RMax, lightDistr.ProfileY.RMax);
                                if isa(registrator, 'PlanarQuadrangle')
                                    registrator.Width = RMax * 2.8;
                                    registrator.Height = RMax * 2.8;
                                end
                        end
                    elseif isa(reqDistr, 'RequiredIntensityDistribution')
                        lightDistr = reqDistr.LightDistribution;
                        mapOptions = IntensityMapOptions.FitToDistribution(lightDistr);
                        project.IntenMapOptions = mapOptions;
                    else
                        mapOptions = IntensityMapOptions.Default;
                    end                    
                    
                    % check cur lens
                    if ~isempty(obj.lenses{obj.curLense, 3})
                        params = obj.lenses{obj.curLense, 3};
                        if isequal(params.reqDistr, reqDistr) && ...
                                isequal(params.source, source)
                            solidParam = params.lens;
                        else
                            solidParam = obj.lenses{obj.curLense, 2};
                        end
                    else
                        solidParam = obj.lenses{obj.curLense, 2};
                    end
                    ch = func(...
                                'Position', getappdata(parentHandle, 'childPosition'),...
                                'frmParent', parentHandle,...
                                'solid', solidParam,...
                                'reqDistr', reqDistr,...
                                'source', source,...
                                'raytrOpt', raytrOpt,...
                                'mapOptions', mapOptions,...
                                'mode', 'Project creation',...
                                'properties', properties);
            end
            uiwait(ch);
            
            % get results and process them
            status = getappdata(parentHandle, 'childResult');
            objRes = getappdata(parentHandle, 'childObject');
                        
            switch step
                case 1
                    obj.addSource(objRes);
                    [isSource, pos] = obj.isSuchSourceIn(objRes);
                    if isSource
                        obj.curSource = pos;
                    else
                        switch objRes
                            case 'StandardSource'
                                newSource = StandardSource.Default;
                            case 'RayFileSource'
                                newSource = RayFileSource.Default;
                        end
                        obj.sources = [obj.sources ; {newSource []}];
                        obj.curSource = size(obj.sources, 1);
                    end
                case 2
                    params = objRes;
                    obj.sources{obj.curSource, 2} = params;
                case 3
                    obj.addReqDistr(objRes);
                case 4
                    params = objRes;
                    obj.requiredDistributions{obj.curReqDistr, 2} = params;
                case 5
                    if ~isempty(objRes)
                        [isLense, pos] = obj.isSuchLenseIn(objRes);
                        if isLense
                            obj.curLense = pos;
                        else
                            lensFunc = str2func([objRes '.DefaultZA']);
                            obj.lenses = [obj.lenses; {objRes lensFunc() []}];
                            obj.curLense = size(obj.lenses, 1);
                        end
                    end
                case 6
                    if isfield(objRes, 'lens')
                        objRes.reqDistr = obj.requiredDistributions{obj.curReqDistr, 1}.clone();
                        if isa(objRes.reqDistr, 'RequiredIlluminanceDistribution')
                            objRes.reqDistr.Registrator = objRes.reqDistr.Registrator.clone();
                        end
                        objRes.source = obj.sources{obj.curSource, 1}.clone();
                        obj.lenses{obj.curLense, 2} = objRes.lens.ZA;
                        obj.lenses{obj.curLense, 3} = objRes;
                    else
                        if ~isequal(obj.lenses{obj.curLense, 2}, objRes)
                            obj.lenses{obj.curLense, 2} = objRes;
                            obj.lenses{obj.curLense, 3} = [];
                        end
                    end
            end            
        end
        
    end
end