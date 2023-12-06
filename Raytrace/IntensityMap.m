classdef (Sealed) IntensityMap < ICloneable & ISerializable
    properties (Access = {?Serializator})
        rays;
        options;
        distr;
        phi;
        psi;
        symDistr;
        rotSymProf;
        interpFunc;
    end
    
    properties (Dependent)
        Rays;
        Options;
        Distr;
        Phi;
        Psi;
        TriPos;
        RotationalSymmetryProfile;
    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods 
        
        function value = get.Rays(obj)
            value = obj.rays;
        end
        
        function value = get.Options(obj)
            value = obj.options.clone();
        end
        
        function value = get.Distr(obj)
            value = obj.symDistr;
        end
        
        function value = get.Phi(obj)
            value = obj.phi;
        end
        
        function value = get.Psi(obj)
            value = obj.psi;
        end
        
        function value = get.TriPos(obj)
            % unpack
            resPhi = obj.phi;
            resPsi = obj.psi;
            
            % make triangles
            i = find( diff(resPsi)~= 0 );
            i = [1; i+1];
            numRings = length(i) - 1;
            i(end+1) = length(resPhi) + 1;
            numPointsInLine = diff(i);
            
            % allocate mem for triangles
            numTriInRings = numPointsInLine(1:end-1) + numPointsInLine(2:end) - 2;
            numTri = sum( numTriInRings );
            triPos = zeros( numTri, 3 );
            triPosIndices = [ 1; cumsum(numTriInRings) + 1 ];
            
            for j = 1 : numRings
                numPoints1 = numPointsInLine(j);
                numPoints2 = numPointsInLine(j+1);
                
                % allocate memory for pointsInds
                numTriInRing = numTriInRings(j);
                pointsInds = zeros(3*numTriInRing, 1);
                
                % go
                pointsInds = twoMeasPointSets2TrianglesMEX( pointsInds, resPhi(i(j) : i(j+1)-1), resPhi(i(j+1) : i(j+2)-1), numPoints1, numPoints2 );
                triPos( triPosIndices(j) : triPosIndices(j+1) - 1, :) = ( reshape(pointsInds, 3, numTriInRing) )' + i(j) - 1;
            end
            
            value = triPos;
        end
        
        function value = get.RotationalSymmetryProfile(obj)
            value = obj.rotSymProf;
        end
       
        function set.Options(obj, value)
            ValueChecker.checkClass(value, 'IntensityMapOptions');
            ValueChecker.checkAndUpdateSize(value, [1 1]);
            
            if obj.options.IsSmoothing ~= value.IsSmoothing ||...
                obj.options.Sig ~= value.Sig || obj.options.NumPoints ~= value.NumPoints ||...
                    obj.options.PsiMax ~= value.PsiMax
                [obj.phi, obj.psi, obj.distr] = DistributionMapProcessor.makeIntensityMap(obj.rays, value);
                [obj.symDistr, obj.rotSymProf] = IntensityMap.applySymmetry(obj.phi, obj.psi, obj.distr, value.Symmetry);
                obj.interpFunc = TriScatteredInterp(obj.phi, obj.psi, obj.symDistr);
            elseif obj.options.Symmetry ~= value.Symmetry
                [obj.symDistr, obj.rotSymProf] = IntensityMap.applySymmetry(obj.phi, obj.psi, obj.distr, value.Symmetry);
                obj.interpFunc = TriScatteredInterp(obj.phi, obj.psi, obj.symDistr);
            end
            obj.options = value;
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PUBLIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access = public)        
        % constructor
        function obj = IntensityMap( varArg1, options, phi, psi, distr )
			% varArg1 - rays
			if nargin == 1 && Serializator.isSerializedObj(varArg1)
				ValueChecker.checkDeserializationPossibility(obj, varArg1);
				obj = Serializator.load(obj, varArg1);
            else
                if(~isempty(varArg1))
                    ValueChecker.checkNumerals(varArg1);
                    varArg1 = ValueChecker.checkAndUpdateSize(varArg1, [inf 7]);          
                    ValueChecker.checkClass(options, 'IntensityMapOptions');
                    ValueChecker.checkAndUpdateSize(options, [1 1]);
                end
                up = options.Up;
                normal = options.Normal;
                lcs.AxisY = up;
                lcs.AxisZ = normal;
                lcs.AxisX = cross(up, normal);
                lcs.Center = [0 0 0];
                
                varArg1 = GeometryConvertor.getLocalRay(varArg1, lcs);
                if nargin == 2
                    [phi, psi, distr] = DistributionMapProcessor.makeIntensityMap(varArg1, options);
                else
                    ValueChecker.checkNumerals(phi);
                    phi = ValueChecker.checkAndUpdateSize(phi, [inf 1]);
                    ValueChecker.checkConstraint(phi, '<=', 2*pi);
                    ValueChecker.checkConstraint(phi, '>=', 0);
                    
                    ValueChecker.checkNumerals(psi);
                    psi = ValueChecker.checkAndUpdateSize(psi, [inf 1]);
                    ValueChecker.checkConstraint(psi, '<=', pi);
                    ValueChecker.checkConstraint(psi, '>=', 0);
                    
                    ValueChecker.checkNumerals(distr);
                    distr = ValueChecker.checkAndUpdateSize(distr, [inf 1]);
                    ValueChecker.checkConstraint(distr, '>=', 0);
                    
                    ValueChecker.checkConstraint(length(distr), '==', length(phi));
                    ValueChecker.checkConstraint(length(distr), '==', length(psi));
                end
                obj.rays = varArg1;
				obj.options = options;
				obj.phi = phi;
				obj.psi = psi;
				obj.distr = distr;
				[obj.symDistr, obj.rotSymProf] = IntensityMap.applySymmetry(phi, psi, distr, options.Symmetry);                
			end % if
            obj.interpFunc = TriScatteredInterp(obj.phi, obj.psi, obj.symDistr);
        end 
        
        function [psi, inten] = getIntensityProfile(obj, phiCur, numPoints)
            if nargin > 2
                ValueChecker.checkNumerals(numPoints);
                ValueChecker.checkAndUpdateSize(numPoints, [1 1]);
                ValueChecker.checkConstraint(phiCur, '>=', 0);                
            else
                numPoints = obj.options.NumPoints;
            end
            
            ValueChecker.checkNumerals(phiCur);
            ValueChecker.checkAndUpdateSize(phiCur, [1 1]);
            ValueChecker.checkConstraint(phiCur, '<=', 2*pi);
            ValueChecker.checkConstraint(phiCur, '>=', 0);
            
            if obj.options.Symmetry == Constants.IntensityMapSymmetries.ROTATIONAL
                F = obj.rotSymProf;
            else
                F = obj.interpFunc;
            end
            psiMin = min(obj.psi);
            psiMax = max(obj.psi);
            dpsi = ( psiMax - psiMin ) / numPoints;
            psi = psiMin : dpsi : psiMax;
            inten = F( repmat(phiCur, size(psi)), psi );
        end
        
        function [inten] = getIntensity(obj, phi, psi)
            ValueChecker.checkNumerals(phi);
            ValueChecker.checkAndUpdateSize(phi, [1 1]);
            ValueChecker.checkConstraint(phi, '<=', 2*pi);
            ValueChecker.checkConstraint(phi, '>=', 0);
            
            ValueChecker.checkNumerals(psi);
            ValueChecker.checkAndUpdateSize(psi, [1 1]);
            ValueChecker.checkConstraint(psi, '<=', pi);
            ValueChecker.checkConstraint(phi, '>=', -pi);
            
            if obj.options.Symmetry == Constants.IntensityMapSymmetries.ROTATIONAL
                F = obj.rotSymProf;
            else
                F = obj.interpFunc;
            end
            inten = F( phi, psi );           
        end
        
        function res = getRRMSE(obj, reqIntenFunc)
            psiMin = reqIntenFunc.PsiMin;
            psiMax = reqIntenFunc.PsiMax;            
            inten = obj.symDistr;
            intenAverLevel = obj.getAverageValue(psiMin, psiMax);
            intenReq = reqIntenFunc.eval(obj.phi, obj.psi);
            
            indGood = obj.psi <= psiMax & obj.psi >= psiMin;
            
            res = sqrt(sum(((inten(indGood) - intenReq(indGood)) .^2 )) / nnz(indGood)) / intenAverLevel * 100;
        end
        
        function res = getMaxValue(obj, psiMin, psiMax)
            if nargin == 1
                psiMin = 0;
                psiMax = pi;
            end            
            inten = obj.symDistr;            
            indGood = obj.psi <= psiMax & obj.psi >= psiMin;
            
            res = max(inten(indGood));
        end
        
        function res = getMinValue(obj, psiMin, psiMax)
            if nargin == 1
                psiMin = 0;
                psiMax = pi;
            end
            inten = obj.symDistr;            
            indGood = obj.psi <= psiMax & obj.psi >= psiMin;
            
            res = min(nonzeros(inten(indGood)));
        end
        
        function res = getAverageValue(obj, psiMin, psiMax)
            if nargin == 1
                psiMin = 0;
                psiMax = pi;
            end            
            inten = obj.symDistr;            
            indGood = obj.psi <= psiMax & obj.psi >= psiMin;
            
            res = sum( inten(indGood) ) / sum(sum(indGood));
        end        
        
        function objNew = clone(obj)
            objNew = IntensityMap(obj.rays, obj.options.clone(), obj.phi, obj.psi, obj.distr);
        end % function
    end % methods
        
    methods (Static = true, Access = private)
        function [symDistr, rotSymFunc] = applySymmetry(phi, psi, distr, symmetry)            
            switch symmetry
                case Constants.IntensityMapSymmetries.QUADRANT                                       
                    distr = IntensityMap.makeUpDownSymmetry(phi, psi, distr);
                    symDistr = IntensityMap.makeLeftRightSymmetry(phi, psi, distr); 
                    rotSymFunc = [];
                case Constants.IntensityMapSymmetries.UD                    
                    symDistr = IntensityMap.makeUpDownSymmetry(phi, psi, distr);
                    rotSymFunc = [];
                case Constants.IntensityMapSymmetries.LR
                    symDistr = IntensityMap.makeLeftRightSymmetry(phi, psi, distr);
                    rotSymFunc = [];
                case Constants.IntensityMapSymmetries.ROTATIONAL
                    % find different psi
                    difPsi = diff(psi);
                    numPoints = nnz(difPsi) + 1;
                    newPsi = zeros(numPoints, 1);
                    newDistr = zeros(numPoints, 1);
                    newPsi(2:end) = psi( find(difPsi ~= 0) + 1 );
                    
                    for i = 1:length(newPsi)
                        n = find( psi == newPsi(i) );
                        newDistr(i) = sum( distr(n) ) / length(n);
                    end
                    
                    rotSymFunc = @(phi, psi) interp1(newPsi, newDistr, psi);
                    symDistr = rotSymFunc([], psi);                    
                otherwise
                    symDistr = distr;
                    rotSymFunc = [];
            end
        end
        
        function distr = makeLeftRightSymmetry(phi, psi, distr)
            F = TriScatteredInterp(phi, psi, distr);
            phiNew = phi;
            m1 = phiNew <= pi & phiNew >= 0;
            m2 = phiNew < 2*pi & phiNew > pi;
            phiNew(m1) = pi - phiNew(m1);
            phiNew(m2) = 3*pi - phiNew(m2);
            distrNew = F(phiNew, psi);
            distr = (distr + distrNew)/2;
        end
        
        function distr = makeUpDownSymmetry(phi, psi, distr)
            F = TriScatteredInterp(phi, psi, distr);
            phiNew = 2*pi - phi;
            distrNew = F(phiNew, psi);
            distr = (distr + distrNew)/2;
        end
    end
end