classdef FreeformBicubSplineSpher < FreeformBicubSplineSurface
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PUBLIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    methods (Access = public)
        
        % constructor
        function obj = FreeformBicubSplineSpher(varArg1, varArg2, iMax, isMirrored, optMatr, varArg3, cs)
            % varArg1 - spline
			if nargin == 1 && Serializator.isSerializedObj(varArg1)
				ValueChecker.checkDeserializationPossibility(obj, varArg1);
				obj = Serializator.load(obj, varArg1);
                obj.optMatr = obj.optMatr .* obj.getFullOptMatr();
            else
                ValueChecker.checkClass(varArg1, 'BicubSpline');
				ValueChecker.checkAndUpdateSize(varArg1, [1 1]);  
				obj.surface = varArg1;   
                
                if nargin == 1
                    iMin = 1;
                    iMax = length(varArg1.X);
                    isMirrored = false;
                    optMatr = [];
                    cs = CoordinateSystem.Global;
                    opticalProperties = OpticalSurfaceProperties.Default;
                elseif nargin == 2
                    iMin = 1;
                    iMax = length(varArg1.X);
                    isMirrored = false;
                    optMatr = [];
                    cs = varArg2;
                    opticalProperties = OpticalSurfaceProperties.Default;
                elseif nargin == 4
                    iMin = varArg2;
                    optMatr = [];
                    cs = CoordinateSystem.Global;
                    opticalProperties = OpticalSurfaceProperties.Default;
                elseif nargin == 5
                    iMin = varArg2;
                    cs = CoordinateSystem.Global;
                    opticalProperties = OpticalSurfaceProperties.Default;
                elseif nargin == 6
                    iMin = varArg2;
                    cs = varArg3;
                    opticalProperties = OpticalSurfaceProperties.Default;                    
                elseif nargin == 7
                    opticalProperties = varArg3;
                    iMin = varArg2;
                else
                    ME = MException('FreeformBicubSplineSpher:FreeformBicubSplineSpher', 'Wrong number of input parameters');
                    throw(ME);
                end
                				
                obj.iMin = iMin;
                obj.iMax = iMax;
                obj.isMirrored = isMirrored;
                obj.cs = cs;
                obj.opticalProperties = opticalProperties;
                if isempty(optMatr)
                    obj.optMatr = obj.getFullOptMatr();
                else
                    obj.optMatr = optMatr .* obj.getFullOptMatr();
                end
			end % if
            if obj.hasCentralPoint()
                obj.updateCentralPoint();
            end
            addlistener(obj.surface, 'SplineChanged', @obj.splineChangedListener);
        end % function
                      
        function bottomSpline = GetBottomContour(obj)
            spline = obj.surface;
            bottomSpline = CubSpline(spline.X, ...
                                     spline.F(end,:) .* sin(spline.Y(end)),...
                                     spline.F_X(end,:) .* sin(spline.Y(end)));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% ICloneable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function objNew = clone(obj)
            objNew = FreeformBicubSplineSpher(obj.surface.clone(), obj.iMin, obj.iMax, obj.isMirrored, obj.optMatr, obj.opticalProperties.clone(), obj.cs.clone());
            objNew.isRegistrator = obj.isRegistrator;
        end % function
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% IDrawable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [pointsCloud, lines] = getPointsCloudAndLines(obj)
            lines = {};
            rspline = obj.surface;
            psi_min = rspline.Y(1);
            psi_max = rspline.Y(end);
%             phi_max = rspline.X(obj.iMax);
            phi_max = rspline.X(end);
            dpsi = (psi_max - psi_min) / 100;
            dphi = phi_max / 100;
            [phi, psi] = meshgrid(0 : dphi : phi_max, psi_min : dpsi : psi_max);
            
            % get r array and calculate x,y,z
            r = rspline.eval(phi, psi);
            x = r .* sin(psi) .* cos(phi);
            y = r .* sin(psi) .* sin(phi);
            z = r .* cos(psi);
            defSize = size(x);
            
            x = reshape(x, [numel(x), 1]);
            y = reshape(y, [numel(y), 1]);
            z = reshape(z, [numel(z), 1]);            
            points = GeometryConvertor.getGlobalPoint([x y z], obj.cs);
            
            pointsCloud(:,:,1) = reshape(points(:,1), defSize);
            pointsCloud(:,:,2) = reshape(points(:,2), defSize);
            pointsCloud(:,:,3) = reshape(points(:,3), defSize);
            
            lines{end+1} = [pointsCloud(end, :, 1)' pointsCloud(end, :, 2)' pointsCloud(end, :, 3)'];
            
            psi = rspline.Y;
%             phi = rspline.X(1:obj.iMax);
            phi = rspline.X;
            % vertical lines
            psi_t = 0 : dpsi: psi_max;
            for j = 1 : length(phi)
                phi_t = repmat(phi(j), size(psi_t));
                r = rspline.eval(phi_t, psi_t);
                x = r .* sin(psi_t) .* cos(phi_t);
                y = r .* sin(psi_t) .* sin(phi_t);
                z = r .* cos(psi_t);
                lines{end+1} = GeometryConvertor.getGlobalPoint([x' y' z'], obj.cs);
            end % for
            
            % horizontal lines
            phi_t = phi(1) : dphi : phi(end);
            for j = 1 : length(psi)
                psi_t = repmat(psi(j), size(phi_t));
                r = rspline.eval(phi_t, psi_t);
                x = r .* sin(psi_t) .* cos(phi_t);
                y = r .* sin(psi_t) .* sin(phi_t);
                z = r .* cos(psi_t);
                lines{end+1} = GeometryConvertor.getGlobalPoint([x' y' z'], obj.cs);
            end % for
        end
        
        function line2D = getProfile2D(obj, phi, numPoints)
            if nargin < 2
                phi = 0;
            end % if
            if nargin < 3
                numPoints = 200;
            end % if
            
            % compute profile
            spline = obj.surface;
            psi = spline.Y(1) : (spline.Y(end) - spline.Y(1)) / (numPoints-1) : spline.Y(end);
            phi = repmat(phi, size(x));
            r = spline.eval(phi, psi);
            x = r .* sin(psi);
            z = r .* cos(psi);
            
            % pack
            line2D = SegLine2D([x' z']);
        end % function

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% IRhinoExportable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function res = getRhinoCommandsString(obj, ~, params)
            npoints = params.NPoints;
            ncurves = params.NCurves;
            spline = obj.surface;
            res = [];
            
            psi_min = spline.Y(1);
            psi_max = spline.Y(end);
            phi_max = spline.X(end);
            dpsi = (psi_max - psi_min) / (npoints - 1);
            if obj.isSurfaceAzimuthalClosed
                dphi = phi_max / ncurves;
                [phi, psi] = meshgrid(0 : dphi : phi_max - dphi, psi_min : dpsi : psi_max);
            else
                dphi = phi_max / (ncurves - 1);
                [phi, psi] = meshgrid(0 : dphi : phi_max, psi_min : dpsi : psi_max);
            end
                        
            % get r array and calculate x,y,z
            r = spline.eval(phi, psi);
            x = r .* sin(psi) .* cos(phi);
            y = r .* sin(psi) .* sin(phi);
            z = r .* cos(psi);
            
            if params.IsMoldable
                xDir = sign(x(2,:) - x(1,:));
                yDir = sign(y(2,:) - y(1,:));
                incAng = deg2rad(2);
                for i = 1 : size(phi,2)
                    for j = 3 : size(psi,1)
                        dx = sign(x(j,i) - x(j-1,i));
                        dy = sign(y(j,i) - y(j-1,i));
                        if (dx ~= xDir(i) || dy ~= yDir(i))
                            lastPoint = [x(j-1,i) y(j-1,i) z(j-1,i)];
                            curPoint = [x(j,i) y(j,i) z(j,i)];
                            curPsi = psi(j, i);
                            curPhi = phi(j, i);
                            
                            lastVec = lastPoint - curPoint;
                            lastR = sqrt(dot(lastVec, lastVec));
                            lastVec = lastVec / lastR;
                            curVec = [sin(curPsi) * cos(curPhi) sin(curPsi) * sin(curPhi) cos(curPsi)];
                            
                            alfa = acos(dot(lastVec, [0 0 1]));
                            betta = acos(dot(lastVec, curVec));
                            
                            incDR = lastR * sin(incAng) * sin(betta) / sin(pi - betta - alfa) / sin(pi - curPsi - incAng);
                            dr = incDR + lastR * sin(alfa) / sin(pi - betta - alfa);
                            newR = r(j,i) + dr;
                            
                            x(j,i) = newR .* sin(curPsi) .* cos(curPhi);
                            y(j,i) = newR .* sin(curPsi) .* sin(curPhi);
                            z(j,i) = newR .* cos(curPsi);
                        end
                    end
                end
            end
            
            for i = 1:ncurves
                %convert points
                point = GeometryConvertor.getGlobalPoint([x(:,i) y(:,i) z(:,i)], obj.cs);
                
                res = [res RhinoExporter.getStringInterpCurve(point(:,1), point(:,2), point(:,3))];
            end
            
            res = [res '_SelCrv\n-_Loft\n_Type\n_Normal\n_Closed=Yes\n_Rebuild\n40\n_Enter\n'];
            res = [res '_SelNone\n_SelCrv\n_Delete\n'];            
        end % function

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% ISTLExportable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function rawTriangs = getTriangles(obj, apprOpt)
            % unpack
            spline = obj.surface;
            maxError = apprOpt.MaxError;
            
            % compute initial points
            numPoints = 100;
            dphi = (spline.X(end) - spline.X(1)) / numPoints;
            dpsi = (spline.Y(end) - spline.Y(1)) / numPoints;
            phi = spline.X(1) : dphi : spline.X(end);
            deltaPhi = spline.X(end) - spline.X(1);
            psi = (spline.Y(1) : dpsi : spline.Y(end))';
            [phi, psi] = meshgrid(phi, psi);
            r = spline.eval(phi, psi);
            x = r .* sin(psi) .* cos(phi);
            y = r .* sin(psi) .* sin(phi);
            z = r .* cos(psi);
            
            % a, b, c
            dlPsiOver1Point =  sqrt((x(2:end, :) - x(1:end-1, :)).^2 + ...
                (y(2:end, :) - y(1:end-1, :)).^2 + ...
                (z(2:end, :) - z(1:end-1, :)).^2);
            dlPsiOver2Points = sqrt((x(3:end, :) - x(1:end-2, :)).^2 + ...
                (y(3:end, :) - y(1:end-2, :)).^2 + ...
                (z(3:end, :) - z(1:end-2, :)).^2);
            dlPhiOver1Point =  sqrt((x(:, 2:end) - x(:, 1:end-1)).^2 + ...
                (y(:, 2:end) - y(:, 1:end-1)).^2 + ...
                (z(:, 2:end) - z(:, 1:end-1)).^2);
            dlPhiOver2Points = sqrt((x(:, 3:end) - x(:, 1:end-2)).^2 + ...
                (y(:, 3:end) - y(:, 1:end-2)).^2 + ...
                (z(:, 3:end) - z(:, 1:end-2)).^2);
            % psi
            aPsi = dlPsiOver1Point(1:end-1,:);
            bPsi = dlPsiOver1Point(2:end,:);
            cPsi = dlPsiOver2Points;
            pPsi = (aPsi + bPsi + cPsi)/2;
            % phi
            aPhi = dlPhiOver1Point(:,1:end-1);
            bPhi = dlPhiOver1Point(:,2:end);
            cPhi = dlPhiOver2Points;
            pPhi = (aPhi + bPhi + cPhi)/2;
            
            % R
            % psi
            RPsi = aPsi .* bPsi .* cPsi ./ ...
                4 ./ sqrt((pPsi .* (pPsi - aPsi) .* (pPsi - bPsi) .* (pPsi - cPsi)));
            RPsi = min(RPsi, [], 2);
            RPsi = [RPsi(1); RPsi; RPsi(end)];
            % phi
            RPhi = aPhi .* bPhi .* cPhi ./ ...
                4 ./ sqrt((pPhi .* (pPhi - aPhi) .* (pPhi - bPhi) .* (pPhi - cPhi)));
            RPhi = min(RPhi, [], 2);
            RPhi(1) = RPhi(2); % suxx
            % deltaLMax by dpsi
            % psi
            dlPsiMax = max(dlPsiOver1Point, [], 2);
            dlPsiMax = [dlPsiMax; dlPsiMax(end)];
            dldpsiPsiMax = dlPsiMax/(psi(2,1) - psi(1,1));
            dpsiErr = sqrt(8 * RPsi * maxError) ./ dldpsiPsiMax;
            % phi
            dlPhiMax = max(dlPhiOver1Point, [], 2);
            dldphiPsiMax = dlPhiMax/(phi(1,2) - phi(1,1));
            dldphiPsiMax(1) = dldphiPsiMax(2); % suxx
            dphiErr = sqrt(8 * RPhi * maxError) ./ dldphiPsiMax;
            
            % make psiNew vector
            psi = psi(:,1); % roll up matrix to vector
            % allocate memory for psiNew vector
            dpsiMin = min(dpsiErr);
            psiNew = zeros(ceil((psi(end) - psi(1)) / dpsiMin) + 1, 1);
            % build psiNew
            psiCur = psi(1);
            psiNew(1) = psiCur;
            curInd = 1;
            dpsiCur = interp1(psi, dpsiErr, psiCur, 'pchip');
            while psiCur < psi(end)
                dpsiNext = interp1(psi, dpsiErr, psiCur + dpsiCur, 'pchip');
                while dpsiNext < dpsiCur
                    dpsiCur = dpsiNext;
                    dpsiNext = interp1(psi, dpsiErr, psiCur + dpsiCur, 'pchip');
                end % while
                psiCur = psiCur + dpsiCur;
                curInd = curInd + 1;
                psiNew(curInd) = psiCur;
                dpsiCur = dpsiNext;
            end % while
            psiNew = psiNew(1:curInd);
            psiNew(end) = psi(end);
            
            % define maximal number of triangles and allocate memory
            dPhiMin = min(sqrt(8 * RPhi * maxError) ./ dldphiPsiMax);
            maxTriangNum = ceil(length(psiNew) * deltaPhi/dPhiMin * 2);
            rawTriangs = zeros(maxTriangNum, 12);
            
            % build triangles
            curTriangNum = 1;
            [points1, phi1] = FreeformBicubSplineSpher.bicubSplineConSec2points(spline, dphiErr(1), psiNew(1));
            for i = 2 : length(psiNew)
                
                % make new line of points
                dphiErrCur = interp1(psi, dphiErr, psiNew(i), 'pchip');
                [points2, phi2] = FreeformBicubSplineSpher.bicubSplineConSec2points(spline, dphiErrCur, psiNew(i));
                
                % pack points
                trianglesI = TriangleSet.twoMeasPointSets2RawTriangles(points1, phi1, points2, phi2);
                numTriangsInRing = size(trianglesI, 1);
                rawTriangs(curTriangNum : numTriangsInRing + curTriangNum - 1, 1:9) = trianglesI;
                curTriangNum = curTriangNum + numTriangsInRing;
                
                % set new line of points to old
                points1 = points2;
                phi1 = phi2;
                
            end % for
            
            % free unused memory
            rawTriangs = rawTriangs(1 : curTriangNum - 1, :);
        end % function
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% I3DEditable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function nodeVector = getNodes(obj)
            phi = obj.surface.X;
            psi = obj.surface.Y;
            R = obj.surface.F;
            if obj.isSurfaceAzimuthalClosed
                numProfiles = length(phi) - 1;
            else
                numProfiles = length(phi);
            end
            numPoints = length(psi);
            nodeVector = zeros(numPoints, numProfiles, 3);
                 
            for i = 1 : numProfiles
                for j = 1 : numPoints
                    curPoint = [ R(j, i) * sin(psi(j)) * cos(phi(i)) ...
                                 R(j, i) * sin(psi(j)) * sin(phi(i)) ...
                                 R(j, i) * cos(psi(j)) ];
                    nodeVector(j, i, :) = GeometryConvertor.getGlobalPoint(curPoint, obj.cs);
                end
            end
        end        
                
    end % methods
    
    methods(Static = true, Access = public)
        
        function str = getCoordinatesType()
            str = 'spher';
        end % function
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PRIVATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access = protected)
        
        %
        function updateCentralPoint(obj, spline)
            if nargin == 1
                spline = obj.surface;
            end
            
            % get 5 points
            delta = spline.Y(2);
            phi = [0 0     pi/2  pi    3*pi/2]';
            psi = [0 delta delta delta delta]';
            r   = spline.eval(phi, psi);
            p   = [r.*sin(psi).*cos(phi)  r.*sin(psi).*sin(phi)  r.*cos(psi)];
            % get tangent vectors and normal
            t1  = p(2, :) - p(4, :);
            t2  = p(3, :) - p(5, :);
            n   = cross(t1, t2);
            % get A, B, C, D
            A = n(1); B = n(2); C = n(3);
            D = -(A * p(1, 1) + B * p(1, 2) + C * p(1, 3));
            
            % get and update unique part
            [f, f_x, f_y, f_xy, x] = spline.getUniquePart(obj.iMin, obj.iMax);
            if length(x) > 1
                f(1, :) = f(1, 1);
%                 f_x(1, :) = 0;
%                 f_y(1, :) = 0;
%                 f_xy(1, :) = 0;
%                 f_y(1, :) = D/C/C * (A * cos(x) + B * sin(x));
%                 f_xy(1, :) = D/C/C * (B * cos(x) - A * sin(x));
            else
                f(1, :) = f(1, 1);
%                 f_x(1, :) = 0;
%                 f_y(1, :) = 0;
%                 f_xy(1, :) = 0;
            end % if
            obj.surface.updateByUniquePart(obj.iMin, obj.isMirrored, f, f_x, f_y, f_xy);
        end        
        
    end % methods
        
    methods (Static = true, Access = private)
               
        %
        function [F_PSI, F_PHIPSI] = getZeroPsiPlaneDerivs(R, R_psi, phi, phiResult) % R - scalar, R_psi and phi are the same size    
            alfa1 = atan(R_psi(1)/R);
            k1 = [cos(alfa1)*cos(phi(1)) cos(alfa1)*sin(phi(1)) sin(alfa1)];
            alfa2 = atan(R_psi(2)/R);
            k2 = [cos(alfa2)*cos(phi(2)) cos(alfa2)*sin(phi(2)) sin(alfa2)];
            N = cross(k1, k2);
            
            tmp = N(3) * R / N(3)^2;
            F_PSI = - tmp * (N(1)*cos(phiResult) + N(2)*sin(phiResult));
            F_PHIPSI = tmp * (N(2)*cos(phiResult) - N(1)*sin(phiResult));
        end
         
        %
        function [points, phi] = bicubSplineConSec2points(spline, dphi, psi)
            if psi ~= 0
                phiMin = spline.X(1);
                phiMax = spline.X(end);
                deltaPhi = phiMax - phiMin;
                numPhi = ceil(deltaPhi / dphi);
                dphi = deltaPhi / numPhi;
                phi = (phiMin : dphi : phiMax)';
                psi = psi(ones(length(phi), 1));
                r = spline.eval(phi, psi);
                x = r .* sin(psi) .* cos(phi);
                y = r .* sin(psi) .* sin(phi);
                z = r .* cos(psi);
                points = [x y z];
                if phiMin == 0 && phiMax == 2*pi
                    points(end, :) = points(1, :);
                end % if
            else
                z = spline.eval(0, 0);
                points = [0 0 z];
                phi = NaN;
            end % if
            
        end % function
            
    end % methods
    
end % classdef