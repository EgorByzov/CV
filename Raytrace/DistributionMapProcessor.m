classdef (Sealed) DistributionMapProcessor
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PUBLIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Static = true, Access = public)
        
        function [x, y, irrad] = makeIlluminanceMap(traces, exitPlane, options)
            % check input
            if ~isa(options, 'IlluminanceMapOptions')
                ME = MException('DistributionMapProcessor:CheckIlluminanceMapOptions', 'Not a valid options object.');
                throw(ME);
            end
            if ~isa(exitPlane, 'PlanarQuadrangle')
                ME = MException('DistributionMapProcessor:CheckExitPlane', 'Not a valid exit plane object. Can be only PlanarQuadrangle.');
                throw(ME);
            end
            
            % make illuminance map
            % unpack parameters
            exitPlaneBorder = exitPlane.Border; 
            lenX = max(exitPlaneBorder(:,1)) - min(exitPlaneBorder(:,1));
            lenY = max(exitPlaneBorder(:,2)) - min(exitPlaneBorder(:,2));
            numPointsX = options.NumPointsX;
            numPointsY = options.NumPointsY;
            sigX = options.SigX;
            sigY = options.SigY;
            isSmoothing = options.IsSmoothing;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % prepare x-y basic grid
            dx = lenX / numPointsX;
            dy = lenY / numPointsY;
            [x, y] = DistributionMapProcessor.makeIrmapXY( lenX, lenY, numPointsX, numPointsY);
            
            if ~isempty(traces)
                ValueChecker.checkNumerals(traces);
                traces = ValueChecker.checkAndUpdateSize(traces, [inf 3]);
            else
                irrad = zeros(size(x));
                return;
            end
            
            if isSmoothing
                
                irrad = makeIrmapMEX(traces, dx*sigX, dy*sigY, x(1,:), y(:,1)', dx, dy);
                irrad = irrad * 1e6;
                
            else
                % choose valid rays extended Illuminance
                %   get pseudo indices
                i = ceil((traces(:,1) + lenX/2) / dx);
                j = ceil((traces(:,2) + lenY/2) / dy);
                %   left only useful rays with indices form 1 to numPoints*e
                ind = find(i > 0   &   i < numPointsX + 1   &   ...
                    j > 0   &   j < numPointsY + 1);
                traces = traces(ind, :);
                i = i(ind);
                j = j(ind);
                numRays = length(ind);
                
                % fill ext Illuminance grid
                irrad = zeros(numPointsY, numPointsX);
                for k = 1 : numRays
                    irrad(j(k), i(k)) = irrad(j(k), i(k)) + traces(k, 3);
                end % for
                
                % convert flux to Illuminance and consider the symmetry
                irrad = irrad / dx / dy * 1e6;
                
            end
        end
        
        function [resPhi, resPsi, inten] = makeIntensityMap(rays, options)
            % check input
            if ~isa(options, 'IntensityMapOptions')
                ME = MException('DistributionMapProcessor:CheckIntensityMapOptions', 'Not a valid options object.');
                throw(ME);
            end
            ValueChecker.checkNumerals(rays);
            rays = ValueChecker.checkAndUpdateSize(rays, [inf 7]);
            
            % make intensity map
            % unpack parameters
            psiMax = options.PsiMax;
            dAlfa = psiMax / options.NumPoints;
            sigma = sin(options.Sig/180*pi * 3) / 3;
            isSmoothing = options.IsSmoothing;
            
            % check sigma
            if sigma > 1/3
                sigma = 1/3;
            elseif sigma < 1e-6
                sigma = 1e-6;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % make phi & psi map
            numRings = ceil(psiMax / dAlfa);
            psi = 0 : psiMax / numRings : psiMax;
            numPointsInLines = ones(length(psi), 1);
            phi = cell(1, numRings + 1);
            phi{1} = 0;
            cosDAlfa = cos(dAlfa);
            for i = 2 : numRings + 1
                sinPsi = sin(psi(i));
                cosPsi = cos(psi(i));
                dphi = acos( (cosDAlfa - cosPsi^2) / sinPsi^2 );
                curNum = ceil( 2*pi/ dphi );
                dphi = 2*pi / curNum;
                phi{i} = 0: dphi : 2*pi - dphi;
                numPointsInLines(i) = length(phi{i});
            end
            
            % prepare params for res(with 2pi)
            inten = zeros(sum(numPointsInLines) + numRings + 1, 1);
            resPhi = zeros(size(inten));
            resPsi = zeros(size(inten));
            % for calculateing( without 2pi )
            pointsPhi = (cell2mat(phi))';
            pointsPsi = zeros(size(pointsPhi));
            solidAngles = zeros(size(pointsPhi)); % for non smoothing
            curPos = 1;
            for i = 1:length(phi)
                pointsPsi(curPos:curPos+numPointsInLines(i)-1) = psi(i);
                if ~isSmoothing && curPos == 1
                    solidAngles(curPos:curPos+numPointsInLines(i)-1) = 2*pi*(1 - cos(psi(i) + dAlfa/2))/ numPointsInLines(i);
                elseif ~isSmoothing
                    solidAngles(curPos:curPos+numPointsInLines(i)-1) = 2*pi*( -cos(psi(i) + dAlfa/2) + cos(psi(i-1) + dAlfa/2) )/ (numPointsInLines(i) + 1);
                end
                
                curPos = curPos+numPointsInLines(i);
            end
            % result points vectors
            P = [ sin(pointsPsi) .* cos(pointsPhi) ...
                sin(pointsPsi) .* sin(pointsPhi) ...
                cos(pointsPsi) ];
            
            if ~isempty(rays)
                ValueChecker.checkNumerals(rays);
                rays = ValueChecker.checkAndUpdateSize(rays, [inf 7]);
                
                if isSmoothing
                    pointsInten = makeIntenMapMEX(rays(:, 4:7), P, sigma, isSmoothing);
                else % no Smoothing
                    cosDAlfa = cos(0.68 * dAlfa);
                    pointsInten = makeIntenMapMEX(rays(:, 4:7), P, cosDAlfa, isSmoothing);
                    pointsInten = pointsInten ./ solidAngles;
                end
            else
                pointsInten = zeros(size(P,1));
            end
            
            
            % duplicate for 2*pi values
            linesPositions = cumsum(numPointsInLines);
            curResPos = 1;
            curPointsPos = 1;
            for i = 1:length(phi)
                nextResPos = curResPos + numPointsInLines(i);
                nextPointsPos = linesPositions(i);
                resPhi(curResPos:nextResPos) = [pointsPhi(curPointsPos:nextPointsPos); 2*pi];
                resPsi(curResPos:nextResPos) = psi(i);
                inten(curResPos:nextResPos) = [pointsInten(curPointsPos:nextPointsPos); pointsInten(curPointsPos)];
                curResPos = nextResPos + 1;
                curPointsPos = nextPointsPos + 1;
            end        
        end
        
    end % methods
    
    
    methods (Static = true, Access = private)
        
        function [x, y] = makeIrmapXY( lenX, lenY, numPointsX, numPointsY)
            dx = lenX / numPointsX;
            dy = lenY / numPointsY;
            [x, y] = meshgrid(-lenX/2 + dx/2 : dx : lenX/2 - dx/2, ...
                -lenY/2 + dy/2 : dy : lenY/2 - dy/2);
        end

    end
    
end % classdef

