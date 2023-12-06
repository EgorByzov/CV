classdef (Sealed) AnalyticalEngine
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PUBLIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Static = true, Access = public)
        
        %%%%%%%%%%%%%%%%%%%%%% PROFILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [psi, r] = computeProfInten2Illum(intenProf, illumProf, R0, nRef, z0)
            effBord = [0.5 1];
            numPoints = 100;
            normConst = AnalyticalEngine.computeNormConstInten2Illum(intenProf, illumProf);
            eff = fzero(@AnalyticalEngine.getProfileInten2Illum, effBord, [], intenProf, illumProf, R0, nRef, z0, normConst, numPoints);
            [~, psi, lp] = AnalyticalEngine.getProfileInten2Illum(eff, intenProf, illumProf, R0, nRef, z0, normConst, numPoints);
            r = lp(:,1);
        end % function
        
        function [psi, r] = computeProfInten2Inten(intenProfIn, intenProfOut, R0, nRef)
            effBord = [0.5 1];
            numPoints = 300;
            normConst = AnalyticalEngine.computeNormConstInten2Inten(intenProfIn, intenProfOut);
            eff = fzero(@AnalyticalEngine.getProfileInten2Inten, effBord, [], intenProfIn, intenProfOut, R0, nRef, normConst, numPoints);
            [~, psi, r] = AnalyticalEngine.getProfileInten2Inten(eff, intenProfIn, intenProfOut, R0, nRef, normConst, numPoints);
        end % function
        
        function [r, z] = computeProfColIllum2Inten(illumProf, intenProf, illumZ, nRef, isConvex, numPoints)
            [ rIn, psiOut ] = AnalyticalEngine.computeCorrFuncIllum2Inten( illumProf, intenProf, numPoints );
            if isConvex
                psiOut = -psiOut;
            end            
            % get lens profile
            [r, z] = ode45(@AnalyticalEngine.diffeqnProfColIllum2Inten, rIn, illumZ, [], rIn, psiOut, nRef);
        end
        
        function [r, z] = computeProfColIllum2Illum(illumProfIn, illumProfOut, illumZIn, illumZOut, nRef, isConvex, numPoints)
            [ rIn, rOut ] = AnalyticalEngine.computeCorrFuncIllum2Illum( illumProfIn, illumProfOut, numPoints );
            if isConvex
                rOut = -rOut;
            end
            [r, z] = ode45(@AnalyticalEngine.diffeqnProfColIllum2Illum, rIn, illumZIn, [], rIn, rOut, illumZOut, nRef);
        end
        
        function [psiOut, R, gam, l] = compute2ProfInten2Illum(intenProf, illumProf, gam0, z0, R0, r0, nRef, hardOut, numPoints)
            if hardOut == 1
                [r, psiOut] = AnalyticalEngine.computeProfInten2Illum(intenProf, illumProf, r0, nRef, z0);
                R = repmat(R0, size(r));
                l = r - R;
                gam = psiOut;
            else
                % unpack
                Rad1 = illumProf.RMin;
                Rad2 = illumProf.RMax;
                
                % ray-corresponding function
                [psi, u] = AnalyticalEngine.computeCorrFuncInten2Illum(intenProf, illumProf, numPoints);
                
                u(isinf(u)) = NaN;
                [uMaxReal, iMax] = max(abs(u));
                Rad = max(Rad1, Rad2);
                u = interp1(psi(1:iMax), u(1:iMax)*abs(Rad/uMaxReal), psi/psi(end)*psi(iMax));
                
                coeffs = polyfit(psi, u, 13);
                
                % l, gam, R
                [psiOut, l_gam_R] = ode45(@AnalyticalEngine.diffeqn2ProfInten2Illum, psi, [r0-R0 gam0 R0]', [],  1, nRef, z0, coeffs, hardOut);
                
                % unpacking
                l = l_gam_R(:,1);
                gam = l_gam_R(:,2);
                R = l_gam_R(:,3);
                
                % control of profiles intersection
                if find(l <= 0, 1)
                    ME = MException('CalculateZA:ProfilesIntersect', ...
                        'Inner profile intersects outer profile. Decrease height of inner surface or increase height of outer surface');
                    throw(ME);
                end
            end
        end
        
        function [phiOut, R, gam, l] = compute2ProfInten2Inten(intenProfIn, intenProfOut, R0, r0, n1, n2, hardOut, alphaMin, numPoints)
            
            % unpack
            phiMin = intenProfIn.PsiMin;
            phiMax = intenProfIn.PsiMax;
            
            % compute betta by phi
            [phi, betta] = AnalyticalEngine.computeCorrFuncInten2Inten(intenProfIn, intenProfOut, numPoints);            
            % gam by phi
            hardInn = 1 - hardOut; % the more work inner surf
            gam = hardInn * betta + (1-hardInn) * phi;
            % compute normal angle, its correction and correction of gam
            alpha = abs(atan((n2*cos(gam) - n1*cos(phi)) ./ (n2*sin(gam) - n1*sin(phi))));
            isNBad = alpha < alphaMin;
            gam(isNBad) = acos(n1/n2*cos(phi(isNBad)));
            % polynom approximation
            coeffsgam = polyfit(phi, gam, 15);            
            % compute R by phi (inner prof)
            dphi = (phiMax - phiMin) / (numPoints - 1);
            phi = phiMin : dphi : phiMax;
            [phiOut, R] = ode45(@AnalyticalEngine.diffeqn2ProfInten2Inten1, phi, R0, [], n1, n2, coeffsgam);  
            
            coeffsR = polyfit(phiOut, R, 15);
            coeffsgam = polyfit(phiOut, gam, 15);
            coeffsbetta = polyfit(phiOut, betta, 15);
            coeffsRDer = polyder(coeffsR);
            coeffsgamDer = polyder(coeffsgam);            
            % compute l by phi (outer prof)
            dphi = (phiMax - phiMin) / (numPoints - 1);
            phi = phiMin : dphi : phiMax;
            [~, l] = ode45(@AnalyticalEngine.diffeqn2ProfInten2Inten2, phi, r0-R0, [], n1, n2, coeffsbetta, coeffsgam, coeffsR, coeffsRDer, coeffsgamDer);
            
            % control of profiles intersection
            if find(l <= 0, 1)
                ME = MException('CalculateZA:ProfilesIntersect', ...
                    'Inner profile intersects outer profile. Decrease height of inner surface or increase height of outer surface');
                throw(ME);
            end
        end
        
        function [innup, innlat, outlat] = computeTIRProfInten2Inten( intenProfIn, intenProfOut, r0, psiBound, innlatConicAng, skirtWidth, nRef )
            % initialize and unpack params
            numPoints = 200;
            polyOrder = 13;
            if isempty(skirtWidth) || ~isfloat(skirtWidth)
                skirtWidth = 0;
            end  %if
            % INNER-UPPER PROFILE
            [xA, zA] = AnalyticalEngine.computeTIRUppProfInten2Inten(intenProfIn, intenProfOut, r0, psiBound, [], nRef);
            innup = [xA zA]; 
            % INNLAT profile
            % computation
            innlatLength = zA(end) / cos(innlatConicAng);
            dt = innlatLength / (numPoints - 1);
            t = 0 : dt : innlatLength;
            x = xA(end) + (innlatLength - t) * sin(innlatConicAng);
            z = zA(end) - (innlatLength - t) * cos(innlatConicAng);
            innlat = [x' z'];            
            % part B and angle gamma
            % upper point
            curProfIn = IntensityProfile(intenProfIn.Func, psiBound, pi/2);
            [psi1, beta1] = AnalyticalEngine.computeCorrFuncIntenInv2Inten(curProfIn, intenProfOut, numPoints);
            phi = pi/2 - psi1;
            gamma = (asin(sin(phi - innlatConicAng) / nRef)) + innlatConicAng;
            Hc = zA(end) - xA(end)*tan(innlatConicAng - pi/2);
            xB = Hc ./ ( tan(phi) - tan(innlatConicAng - pi/2) );
            zB = xB .* tan(phi);
            rB = sqrt(xB.^2 + zB.^2);            
            % OUTLAT
            % norm constant
            betaC = asin(sin(beta1) / nRef);
            s1 = [-sin(betaC) cos(betaC)];
            % coefficients, part C
            coeffsRB = polyfit (phi, rB, polyOrder);
            coeffsGamma = polyfit(phi, gamma, polyOrder);
            coeffsS1x = polyfit(phi, s1(:, 1), polyOrder);
            coeffsS1z = polyfit(phi, s1(:, 2), polyOrder);
            coeffsRBdPhi = polyder(coeffsRB);
            coeffsGammadPhi = polyder(coeffsGamma);
            [~, l] = ode45(@AnalyticalEngine.diffeqnTIRProfInten2Inten2, phi, skirtWidth / cos(gamma(1)), [], coeffsRB, coeffsGamma,...
                                                            coeffsRBdPhi, coeffsGammadPhi, coeffsS1x, coeffsS1z, nRef);
            xC =  xB + l .* cos(gamma);
            zC =  zB + l .* sin(gamma);
            outlat = [xC zC];
        end
               
        function [xA, zA] = computeTIRUppProfInten2Inten(intenProfIn, intenProfOut, r0, psiBorders, betaBorders, nRef)
            
            % initialize and unpack params
            numPoints = 200;
            if isempty(psiBorders)
                psiMin = intenProfIn.PsiMin;
                psiMax = intenProfIn.PsiMax;
            elseif length(psiBorders) == 1
                psiMin = 0;
                psiMax = psiBorders;
            else
                psiMin = psiBorders(1);
                psiMax = psiBorders(2);
            end % if
            if isempty(betaBorders)
                betaMin = intenProfOut.PsiMin;
                betaMax = intenProfOut.PsiMax;
            elseif length(betaBorders) == 1
                betaMin = 0;
                betaMax = betaBorders;
            else
                betaMin = betaBorders(1);
                betaMax = betaBorders(2);
            end % if
            % INNER-UPPER PROFILE
            % correspondence-function
            curProfIn = IntensityProfile(intenProfIn.Func, psiMin, psiMax);
            curProfOut = IntensityProfile(intenProfOut.Func, betaMin, betaMax);
            [beta0, beta1] = AnalyticalEngine.computeCorrFuncInten2Inten(curProfIn, curProfOut, numPoints);
            betaA = asin(sin(beta1) / nRef);
            % INNUP
            [psi, rA] = ode45(@AnalyticalEngine.diffeqnTIRProfInten2Inten1, beta0, r0, [], beta0, betaA, nRef);
            k = find(isnan(rA), 1, 'first');
            if ~isempty(k)
                psi = psi(1:k-1);
                rA = rA(1:k-1);
            end
            xA = rA .* sin(psi);
            zA = rA .* cos(psi);
        end % function
        
        function [innup, innlat, outlat] = computeTIRProfInten2Illum( intenProf, illumProf, z0, r0, psiBound, innlatConicAng, skirtWidth, nRef )
            % initialize and unpack params
            numPoints = 200;
            polyOrder = 13;
            if isempty(skirtWidth) || ~isfloat(skirtWidth)
                skirtWidth = 0;
            end            
            % INNER-UPPER PROFILE
            % correspondence-function
            curIntenProf = IntensityProfile(intenProf.Func, 0, psiBound);
            [beta0, ro] = AnalyticalEngine.computeCorrFuncInten2Illum(curIntenProf, illumProf, numPoints);
            betaA = asin(sin(atan(ro / z0)) / nRef);            
            % INNUP
            [~, rA] = ode45(@AnalyticalEngine.diffeqnTIRProfInten2Illum1, beta0, r0, [], beta0, betaA, nRef);
            xA = rA .* sin(beta0);
            zA = rA .* cos(beta0);
            innup = [xA zA];            
            % INNLAT profile
            % computation
            innlatLength = zA(end) / cos(innlatConicAng);
            dt = innlatLength / (numPoints - 1);
            t = 0 : dt : innlatLength;
            x = xA(end) + (innlatLength - t) * sin(innlatConicAng);
            z = zA(end) - (innlatLength - t) * cos(innlatConicAng);
            innlat = [x' z'];            
            % part B and angle gamma
            % upper point
            curIntenProf = IntensityProfile(intenProf.Func, psiBound, pi/2);
            [psi1, ro] = AnalyticalEngine.computeCorrFuncIntenInv2Illum(curIntenProf, illumProf, numPoints);
            phi = pi/2 - psi1;
            gamma = (asin(sin(phi - innlatConicAng) / nRef)) + innlatConicAng;
            Hc = zA(end) - xA(end)*tan(innlatConicAng - pi/2);
            xB = Hc ./ ( tan(phi) - tan(innlatConicAng - pi/2) );
            zB = xB .* tan(phi);
            rB = sqrt(xB.^2 + zB.^2);            
            % OUTLAT
            % norm constant
            betaC = asin(sin(atan(ro / z0)) / nRef);
            s1 = [-sin(betaC) cos(betaC)];
            % coefficients, part C
            coeffsRB = polyfit (phi, rB, polyOrder);
            coeffsGamma = polyfit(phi, gamma, polyOrder);
            coeffsS1x = polyfit(phi, s1(:, 1), polyOrder);
            coeffsS1z = polyfit(phi, s1(:, 2), polyOrder);
            coeffsRBdPhi = polyder(coeffsRB);
            coeffsGammadPhi = polyder(coeffsGamma);
            [~, l] = ode45(@AnalyticalEngine.diffeqnTIRProfInten2Illum2, phi, skirtWidth / cos(gamma(1)), [], coeffsRB, coeffsGamma,...
                coeffsRBdPhi, coeffsGammadPhi, coeffsS1x, coeffsS1z, nRef);
            xC =  xB + l .* cos(gamma);
            zC =  zB + l .* sin(gamma);
            outlat = [xC zC];
        end
        
        function [innup, innlat, outlat] = computeTIRCol(r0, psiBound, innlatConicAng, skirtWidth, nRef)
            % initialize and unpack params
            numPoints = 200;
            if isempty(skirtWidth) || ~isfloat(skirtWidth)
                skirtWidth = 0;
            end
            
            % INNER-UPPER PROFILE
            innup = AnalyticalEngine.getColProfileCart(r0, psiBound, nRef, numPoints);
            
            % INNER LATERAL PROFILE
            innlatLength = innup(end,2) / cos(innlatConicAng);
            dt = innlatLength / (numPoints - 1);
            t = 0 : dt : innlatLength;
            x = innup(end,1) + (innlatLength - t) * sin(innlatConicAng);
            z = innup(end,2) - (innlatLength - t) * cos(innlatConicAng);
            innlat = [x' z'];
            
            % TIR PROFILE
            beta = atan(innlat(:,2) ./ innlat(:,1));
            phi = asin(1/nRef * sin(beta - innlatConicAng));
            gamma = innlatConicAng + phi;
            Psi0 = innlat(1, 1) + skirtWidth * nRef / cos(gamma(1));
            H = skirtWidth * tan(gamma(1));
            r = (Psi0 - sqrt(innlat(:,1).^2 + innlat(:,2).^2) - H * nRef + innlat(:,2) * nRef) ./ ...
                nRef ./ (1 - sin(gamma));
            x = innlat(:,1) + r .* cos(gamma);
            z = innlat(:,2) + r .* sin(gamma);
            outlat = [x z];            
        end
        
        function surfs = computeFlatCol2(h0, hRelief, hBase, alphaInc, rRefract, rTIR, rMax, nRef, angleMod)
            if nargin == 8
                angleMod = 0;
            end % if
            
            % compute refractive part
            numPoints = 10;
            % sec 1
            psiMax = acos( (h0 + hRelief) / (h0 + nRef * hRelief) );
            prof = AnalyticalEngine.getColProfileCart(h0, psiMax, nRef, numPoints);
            surfs{1} = AxisymInterpLineCart(prof);
            while true
                % get straight segment
                x1 = surfs{end}.Profile(end, 1);
                z1 = surfs{end}.Profile(end, 2);
                x2 = x1 + (z1 - h0) * tan(alphaInc);
                z2 = h0;
                if x2 > rRefract
                    break;
                end % if
                surfs{end+1} = AxisymSegLineCart(SegLine2D([x1 z1; x2 z2]));
                % get refractive segment
                r = sqrt(x2^2 + z2^2);
                psiMin = atan(x2 / z2);
                psiMax = acos((h0 + hRelief) / (r + nRef * hRelief));
                prof = AnalyticalEngine.getColProfileCart(r, [psiMin psiMax], nRef, numPoints);
                surfs{end+1} = AxisymInterpLineCart(prof);
            end % while
            % compute TIR part
            x1 = surfs{end}.Profile(end, 1);
            z1 = surfs{end}.Profile(end, 2);
            iter = 1;
            numAngleMod = length(angleMod);
            while x1 < rTIR
                x2 = x1 + (z1 - h0) * tan(alphaInc);
                z2 = h0;
                if x2 >= rTIR
                    break;
                end % if
                psi = atan(x2 / z2);
                gamma = asin( cos(psi + alphaInc) / nRef ) + alphaInc;
                bettaInc = pi/4 - gamma/2;
                bettaInc = bettaInc + angleMod(mod(iter, numAngleMod) + 1);    % modificator
                x3 = x2 + hRelief * tan(bettaInc);
                z3 = h0 + hRelief;
                if x3 > rMax
                    break;
                end % if
                surfs{end+1} = AxisymSegLineCart(SegLine2D([x1 z1; x2 z2; x3 z3]));
                x1 = surfs{end}.Profile.Points(end, 1);
                z1 = surfs{end}.Profile.Points(end, 2);
                iter = iter + 1;
            end % while
            % final
            surfs{end+1} = AxisymSegLineCart(SegLine2D(...
                [   x1   z1;
                    rMax h0+hRelief;
                    rMax h0+hRelief+hBase;
                    0    h0+hRelief+hBase...
                                                ]));
        end % function
        
        function [prof1, segProf] = computeFlatCol(r0, hMax, psiBound, alfaMin, RMax, nRef, colAlfa, secondColAlfa, numTriangles)
            if (nargin <= 6)
                numTriangles = 1;
            end
            triangDX = hMax * sin(alfaMin);
            
            % compute first prof
            numPoints = 20;
            psiMax = acos( (r0+hMax)/(r0+nRef*hMax) );
            if psiMax > psiBound
                psiMax = psiBound;
            end % if
            prof1 = AnalyticalEngine.getColProfileCart(r0, psiMax, nRef, numPoints);
            while psiMax < psiBound
                psiMin = atan((prof1(end, 1) + triangDX) / r0);
                if psiMin > psiBound
                    break;
                end % if
                r = sqrt((prof1(end, 1) + triangDX) .^ 2 + r0 .^ 2);
                psiMax = acos((r0 + hMax) / (r + nRef * hMax));
                if psiMax > psiBound
                    psiMax = psiBound;
                end % if
                profAdd = AnalyticalEngine.getColProfileCart(r, [psiMin psiMax], nRef, numPoints);
                prof1 = [prof1; profAdd];
            end % while
            
            % compute second profs triangle part
%             alfa2 = asin( sin(pi/2 - psiBoundCur - alfaMin) / nRef );
%             betta = pi/4 + (alfaMin + alfa2)/2 - colAlfa/2;
%             x1 = hMax * tan(alfaMin);
%             x2 = x1 + hMax / tan(betta);
%             triSegProf = [0     hMax;
%                           x1    0;
%                           x2    hMax];
%             if (nargin > 6)
%                 betta2 = pi/4 + (alfaMin + alfa2)/2 - secondColAlfa/2;
%                 x2_2 = x1 + hMax / tan(betta2);
%                 triSegProf2 = [0     hMax;
%                               x1    0;
%                               x2_2    hMax];
%             end
%                   figure;plot(prof1(:,1),prof1(:,2)); axis equal; ylim([0 20]);    
            % make second profile
            segProf = [];            
            xCur = prof1(end,1);
            iter = 1;
            while xCur < RMax
                curPsi = atan( xCur / (r0 + hMax/2) );
                alfa2 = asin( sin(pi/2 - curPsi - alfaMin) / nRef );
                if numTriangles == 1 || mod(iter,numTriangles) == 1
                    betta = pi/4 + (alfaMin + alfa2)/2 - colAlfa/2;
                    x1 = hMax * tan(alfaMin);
                    x2 = x1 + hMax / tan(betta);
                    segProf = [segProf;
                               xCur     hMax + r0;
                               x1+xCur  r0;
                               x2+xCur  hMax + r0];
                else
                    curTriNum = mod(iter,numTriangles);
                    dAlfa = (secondColAlfa - colAlfa)/numTriangles * (curTriNum-1) + colAlfa;
                    betta = pi/4 + (alfaMin + alfa2)/2 - dAlfa/2;
                    x1 = hMax * tan(alfaMin);
                    x2 = x1 + hMax / tan(betta);
                    segProf = [segProf;
                               xCur     hMax + r0;
                               x1+xCur  r0;
                               x2+xCur  hMax + r0];
                end
                xCur = segProf(end,1);
                iter = iter + 1;
            end                             
        end
        
        function [ r, irrad ] = computeIllumAfterTIR( innup, innlat, outlat, inDistr, radius )
            % unpack params
%             sourceFlux = source.Flux;
            
            % compute
            psiUp = atan(innup(:,1) ./ innup(:,2));
            psiUpMid = (psiUp(2:end) + psiUp(1:end-1)) / 2;
            rUp = innup(:,1);
            rUpMid = interp1(psiUp, rUp, psiUpMid);
            intenUpMid = inDistr.eval(psiUpMid);
            deltaPsiUp = psiUp(2:end) - psiUp(1:end-1);
            deltaRUp = rUp(2:end) - rUp(1:end-1);
            irradUpMid = intenUpMid .* sin(psiUpMid) ./ rUpMid .* abs(deltaPsiUp ./ deltaRUp);
            % lateral section
            psiLat = atan(innlat(:,1) ./ innlat(:,2));
            psiLatMid = (psiLat(2:end) + psiLat(1:end-1)) / 2;
            rLat = outlat(:,1);
            rLatMid = interp1(psiLat, rLat, psiLatMid);
            intenLatMid = inDistr.eval(psiLatMid);
            deltaPsiLat = psiLat(2:end) - psiLat(1:end-1);
            deltaRLat = rLat(2:end) - rLat(1:end-1);
            irradLatMid = intenLatMid .* sin(psiLatMid) ./ rLatMid .* abs(deltaPsiLat ./ deltaRLat);
            % unite r and irrad
            rGap = rLatMid(1) - rUpMid(end);
            r = [0; rUpMid; rUpMid(end) + rGap * 0.01; ...
                rLatMid(1) - rGap * 0.01; rLatMid; radius];
            irrad = [irradUpMid(1); irradUpMid; 0; 0; irradLatMid; irradLatMid(end)];
        end
                
        function [psi, r] = computeReflProfDownInten2Illum(intenProf, illumProf, z0, R0, numPoints)
            
            % unpack
            Rad1 = illumProf.RMin;
            Rad2 = illumProf.RMax;
            
            % ray-corresponding function
            [psi, u] = AnalyticalEngine.computeCorrFuncInten2Illum(intenProf, illumProf, numPoints);
            
            u(isinf(u)) = NaN;
            [uMaxReal, iMax] = max(abs(u));
            Rad = max(Rad1, Rad2);
            u = interp1(psi(1:iMax), u(1:iMax)*abs(Rad/uMaxReal), psi/psi(end)*psi(iMax));
            
            coeffs = polyfit(psi, u, 13);
            
            % psi, r
            [psi, r] = ode45(@AnalyticalEngine.diffeqnReflProfDownInten2Illum, psi, R0, [], z0, coeffs);
            psi = pi - psi;
            
        end % function
       
        function [psi, r] = computeReflProfUpInten2Illum(intenProf, illumProf, z0, R0, alphaMin, alphaMax, numPoints)
                % ray-corresponding function
                [psi, u] = AnalyticalEngine.computeReflCorrFuncUpInten2Illum(intenProf, illumProf, z0, alphaMin, alphaMax, numPoints);
                u = -u; % convex mirror
                
                coeffs = polyfit(psi, u, 13);
                
                % psi, r
                [psi, r] = ode45(@AnalyticalEngine.diffeqnReflProfUpInten2Illum, psi, R0/cos(alphaMin), [], z0, coeffs);                
        end % function

        function [psi, r] = computeReflProfUpInten2Inten(intenProfIn, intenProfOut, R0, alphaMin, alphaMax, numPoints)
            % ray-corresponding function
            [psi, beta] = AnalyticalEngine.computeReflCorrFuncUpInten2Inten(intenProfIn, intenProfOut, alphaMin, alphaMax, numPoints);

            coeffs = polyfit(psi, beta, 13);
            
            % psi, r
            [psi, r] = ode45(@AnalyticalEngine.diffeqnReflProfUpInten2Inten, psi, R0/cos(alphaMin), [], coeffs);
        end % function
        
        function [psi, r] = computeReflProfDownInten2Inten(intenProfIn, intenProfOut, R0, numPoints)
            
            % ray-corresponding function
            [psi, beta] = AnalyticalEngine.computeCorrFuncInten2Inten(intenProfIn, intenProfOut, numPoints);
                        
            coeffs = polyfit(psi, beta, 13);
            
            % psi, r
            [psi, r] = ode45(@AnalyticalEngine.diffeqnReflProfDownInten2Inten, psi, R0, [], coeffs);
            psi = pi - psi;
            
        end % function
        
        function [psi, r] = computeReflCol(R0, alphaMin, alphaMax, numPoints)
            p = R0 * (1 / cos(alphaMin) - 1);
            psi = alphaMin : (alphaMax - alphaMin) / (numPoints - 1) : alphaMax;
            
            r = p ./ ( 1 - cos(psi) );            
        end % function
        
        %%%%%%%%%%%%%%%%%%%%%% NORMALIZATION CONSTANTS %%%%%%%%%%%%%%%%%%%%
        
        function normConst = computeNormConstInten2Illum(intenProf, illumProf)
            sourceFlux = IntensityAxisym(intenProf).computeFlux();
            targetFlux = IlluminanceAxisym(illumProf).computeFlux();
            normConst = sourceFlux / targetFlux;
        end % function
        
        function normConst = computeNormConstInten2Inten(intenProfIn, intenProfOut)
            sourceFlux = IntensityAxisym(intenProfIn).computeFlux();
            targetFlux = IntensityAxisym(intenProfOut).computeFlux();
            normConst = sourceFlux / targetFlux;
        end % function
        
        function normConst = computeNormConstIllum2Inten(illumProf, intenProf)
            sourceFlux = IlluminanceAxisym(illumProf).computeFlux();
            targetFlux = IntensityAxisym(intenProf).computeFlux();
            normConst = sourceFlux / targetFlux;
        end
        
        function normConst = computeNormConstIllum2Illum(illumProfIn, illumProfOut)
            sourceFlux = IlluminanceAxisym(illumProfIn).computeFlux();
            targetFlux = IlluminanceAxisym(illumProfOut).computeFlux();
            normConst = sourceFlux / targetFlux;
        end
        
        %%%%%%%%%%%%%%%%%%%%%% CORRESPONDING FUNCTIONS %%%%%%%%%%%%%%%%%%%%
        
        function [rIn, psiOut] = computeCorrFuncIllum2Inten( illumProf, intenProf, numPoints )
            % unpack
            rIn1 = illumProf.RMin;
            rIn2 = illumProf.RMax;
            psiOut1 = intenProf.PsiMin;
            psiOut2 = intenProf.PsiMax;
            
            % check if out func is bad
            res = fminbnd(@(alfa) intenProf.eval(alfa), psiOut1, psiOut2);
            if intenProf.eval(res) < 1e-6
                ME = MException('CalculateZA:BadEvaluationFunction', ...
                    'Zero or negative intensity function value is not appropriate');
                throw(ME);
            end
            
            % compute norm const
            normConst = AnalyticalEngine.computeNormConstIllum2Inten(illumProf, intenProf);
            % obtain r_out2(r) dependence
            rIn = rIn1 : (rIn2 - rIn1)/ (numPoints - 1) : rIn2;
            [rIn, cosPsiOut] = ode45(@AnalyticalEngine.diffeqnCorrFuncIllum2Inten, rIn, cos(intenProf.PsiMin), [], illumProf, intenProf, normConst);
            psiOut = acos( cosPsiOut );
            % remove bad values
            iMax = find( isinf(psiOut) | isnan(psiOut), 1, 'first' );
            if ~isempty(iMax)
                iMax = iMax - 1;
                psiOut = interp1( (rIn(1:iMax) - rIn(1))/(rIn(iMax) - rIn(1))*(rIn2 - rIn1) + rIn1,...
                    (psiOut(1:iMax) - psiOut(1))/(psiOut(iMax) - psiOut(1))*(psiOut2 - psiOut1) + psiOut1,...
                    (rIn - rIn(1))/(rIn(end) - rIn(1))*(rIn2 - rIn1) + rIn1 );
            end
        end
        
        function [rIn, rOut] = computeCorrFuncIllum2Illum( illumProfIn, illumProfOut, numPoints )
            %unpack
            rIn1 = illumProfIn.RMin;
            rIn2 = illumProfIn.RMax;
            rOut1 = illumProfOut.RMin;
            rOut2 = illumProfOut.RMax;
            
            % check if out func is bad
            res = fminbnd(@(r)illumProfOut.eval(r), rOut1, rOut2);
            if illumProfOut.eval(res) < 1e-6
                ME = MException('CalculateZA:BadEvaluationFunction', ...
                    'Zero or negative irradiance function value is not appropriate');
                throw(ME);
            end   
            
            % compute norm const
            normConst = AnalyticalEngine.computeNormConstIllum2Illum(illumProfIn, illumProfOut);            
            % obtain r_out2(r) dependence
            rIn = rIn1 : (rIn2 - rIn1)/ (numPoints - 1) : rIn2;
            [rIn, rOut] = ode45(@AnalyticalEngine.diffeqnCorrFuncIllum2Illum, rIn, illumProfOut.RMin, [], illumProfIn, illumProfOut, normConst);
            rOut = sqrt(rOut);            
            % remove bad values
            rOut(isinf(rOut)) = NaN;
            [~, iMax] = max(abs(rOut));
            rOut = interp1( (rIn(1:iMax) - rIn(1))/(rIn(iMax) - rIn(1))*(rIn2 - rIn1) + rIn1,...
                            (rOut(1:iMax) - rOut(1))/(rOut(iMax) - rOut(1))*(rOut2 - rOut1) + rOut1,...
                            (rIn - rIn(1))/(rIn(end) - rIn(1))*(rIn2 - rIn1) + rIn1 );
        end
        
        function [psiIn, rOut] = computeCorrFuncInten2Illum( intenProf, illumProf, numPoints )
            % unpack
            psi1 = intenProf.PsiMin;
            psi2 = intenProf.PsiMax;
            rad1 = illumProf.RMin;
            rad2 = illumProf.RMax;
            
            % check if out func is bad
            res = fminbnd(@(r)illumProf.eval(r), rad1, rad2);
            if illumProf.eval(res) < 1e-6
                ME = MException('CalculateZA:BadEvaluationFunction', ...
                    'Zero or negative irradiance function value is not appropriate');
                throw(ME);
            end            
            % compute norm const
            normConst = AnalyticalEngine.computeNormConstInten2Illum(intenProf, illumProf);            
            % compute ray-corresponding function            
            dpsi = (psi2 - psi1) / (numPoints - 1);
            psiIn = psi1 : dpsi : psi2;
            [psiIn, ro2] = ode45(@AnalyticalEngine.diffeqnCorrFuncInten2Illum, psiIn, (illumProf.RMin)^2, [], intenProf, illumProf, normConst);
            rOut = sqrt(ro2);            
            % remove bad values
            rOut(isinf(rOut)) = NaN;
            [~, iMax] = max(abs(rOut));
            rOut = interp1( (psiIn(1:iMax) - psiIn(1))/(psiIn(iMax) - psiIn(1))*(psi2 - psi1) + psi1,...
                (rOut(1:iMax) - rOut(1))/(rOut(iMax) - rOut(1))*(rad2 - rad1) + rad1,...
                (psiIn - psiIn(1))/(psiIn(end) - psiIn(1))*(psi2 - psi1) + psi1);
        end
        
        function [psiIn, rOut] = computeCorrFuncIntenInv2Illum( intenProf, illumProf, numPoints )
            % unpack
            psi1 = intenProf.PsiMin;
            psi2 = intenProf.PsiMax;
            rad1 = illumProf.RMin;
            rad2 = illumProf.RMax;
            
            % check if out func is bad
            res = fminbnd(@(r)illumProf.eval(r), rad1, rad2);
            if illumProf.eval(res) < 1e-6
                ME = MException('CalculateZA:BadEvaluationFunction', ...
                    'Zero or negative irradiance function value is not appropriate');
                throw(ME);
            end            
            % compute norm const
            normConst = AnalyticalEngine.computeNormConstInten2Illum(intenProf, illumProf);            
            % compute ray-corresponding function            
            dpsi = (psi2 - psi1) / (numPoints - 1);
            psiIn = psi2 : -dpsi : psi1;
            [psiIn, ro2] = ode45(@AnalyticalEngine.diffeqnCorrFuncInten2Illum, psiIn, illumProf.RMin, [], intenProf, illumProf, normConst);
            rOut = sqrt(abs(ro2));            
            % remove bad values
            rOut(isinf(rOut)) = NaN;
            [~, iMax] = max(abs(rOut));
            rOut = interp1( (psiIn(1:iMax) - psiIn(1))/(psiIn(iMax) - psiIn(1))*(psi2 - psi1) + psi1,...
                (rOut(1:iMax) - rOut(1))/(rOut(iMax) - rOut(1))*(rad2 - rad1) + rad1,...
                (psiIn - psiIn(1))/(psiIn(end) - psiIn(1))*(psi2 - psi1) + psi1);
        end
        
        function [psiIn, psiOut] = computeCorrFuncInten2Inten(intenProfIn, intenProfOut, numPoints)
            % unpack
            psi1 = intenProfIn.PsiMin;
            psi2 = intenProfIn.PsiMax;
            betta1 = intenProfOut.PsiMin;
            betta2 = intenProfOut.PsiMax;
            
            % check if out func is bad
            res = fminbnd(@(alfa) intenProfOut.eval(alfa), betta1, betta2);
            if intenProfOut.eval(res) < 1e-6
                ME = MException('CalculateZA:BadEvaluationFunction', ...
                    'Zero or negative intensity function value is not appropriate');
                throw(ME);
            end            
            % compute normConst
            normConst = AnalyticalEngine.computeNormConstInten2Inten(intenProfIn, intenProfOut);            
            % compute ray-corresponding function            
            dpsi = (psi2-psi1) / (numPoints - 1);
            psiIn = psi1 : dpsi : psi2;
            [psiIn, lp] = ode45(@AnalyticalEngine.diffeqnCorrFuncInten2Inten, psiIn, cos(intenProfOut.PsiMin), [], intenProfIn, intenProfOut, normConst);
            psiOut = acos(lp);            
            % remove bad values
            iMax = find( isinf(psiOut) | isnan(psiOut), 1, 'first' );
            if ~isempty(iMax)
                iMax = iMax - 1;
                psiOut = interp1( (psiIn(1:iMax) - psiIn(1))/(psiIn(iMax) - psiIn(1))*(psi2 - psi1) + psi1,...
                    (psiOut(1:iMax) - psiOut(1))/(psiOut(iMax) - psiOut(1))*(betta2 - betta1) + betta1,...
                    (psiIn - psiIn(1))/(psiIn(end) - psiIn(1))*(psi2 - psi1) + psi1 );
            end
        end        
        
        function [psiIn, psiOut] = computeCorrFuncIntenInv2Inten(intenProfIn, intenProfOut, numPoints)
            % unpack
            psi1 = intenProfIn.PsiMin;
            psi2 = intenProfIn.PsiMax;
            betta1 = intenProfOut.PsiMin;
            betta2 = intenProfOut.PsiMax;
            
            % check if out func is bad
            res = fminbnd(@(alfa) intenProfOut.eval(alfa), betta1, betta2);
            if intenProfOut.eval(res) < 1e-6
                ME = MException('CalculateZA:BadEvaluationFunction', ...
                    'Zero or negative intensity function value is not appropriate');
                throw(ME);
            end            
            % compute normConst
            normConst = -AnalyticalEngine.computeNormConstInten2Inten(intenProfIn, intenProfOut);            
            % compute ray-corresponding function            
            dpsi = (psi2-psi1) / (numPoints - 1);
            psiIn = psi2 : -dpsi : psi1;
            [psiIn, lp] = ode45(@AnalyticalEngine.diffeqnCorrFuncInten2Inten, psiIn, cos(intenProfOut.PsiMin), [], intenProfIn, intenProfOut, normConst);
            psiOut = acos(abs(lp));            
            % remove bad values
            iMax = find( isinf(psiOut) | isnan(psiOut), 1, 'first' );
            if ~isempty(iMax)
                iMax = iMax - 1;
                psiOut = interp1( (psiIn(1:iMax) - psiIn(1))/(psiIn(iMax) - psiIn(1))*(psi2 - psi1) + psi1,...
                    (psiOut(1:iMax) - psiOut(1))/(psiOut(iMax) - psiOut(1))*(betta2 - betta1) + betta1,...
                    (psiIn - psiIn(1))/(psiIn(end) - psiIn(1))*(psi2 - psi1) + psi1 );
            end
        end
        
        function [psi, rOut] = computeReflCorrFuncUpInten2Illum( intenProf, illumProf, z0, alphaMin, alphaMax, numPoints )
            
            % unpack
            psi1 = intenProf.PsiMin;
            rad1 = illumProf.RMin;
            rad2 = illumProf.RMax;
            
            % check if out func is bad
            res = fminbnd(@(r)illumProf.eval(r), rad1, rad2);
            if illumProf.eval(res) < 1e-6
                ME = MException('CalculateZA:BadEvaluationFunction', ...
                    'Zero or negative irradiance function value is not appropriate');
                throw(ME);
            end % if
            
            % preparing reqIllumNew
            dpsiMir = (alphaMax - alphaMin) / (numPoints - 1);
            psiMir = alphaMin : dpsiMir : alphaMax; 
            dx = (rad2 - rad1) / (numPoints - 1);
            x = rad1 : dx : rad2;
            dpsiDir = (alphaMin - psi1) / numPoints;
            psiDir = psi1 + dpsiDir/2 : dpsiDir : alphaMin - dpsiDir/2;
            xDir = z0 * tan(psiDir);
            illumDir = IntensityAxisym(intenProf).scaleEval([], psiDir, 1) .* (cos(psiDir)).^3 / z0^2 * 1e6;
            xDir = [z0 * tan(psi1) xDir z0 * tan(alphaMin)];
            illumDir = [illumDir(1) illumDir illumDir(end)]; 
            illumDirInterp = Interp2D(xDir, illumDir, 'cubic');
            illumDir = illumDirInterp.eval(x);
            illumVal = IlluminanceAxisym(illumProf).scaleEval(0, x, 1);
            illumReqNew = illumVal - illumDir;
            illumReqNew(illumReqNew < 0) = 0;
            illumReqNew(isnan(illumReqNew)) = 0;
            firstNonZero = find(illumReqNew ~= 0, 1, 'first');
            lastNonZero = find(illumReqNew ~= 0, 1, 'last');
            illumReqNew = illumReqNew(firstNonZero : lastNonZero);
            x = x(firstNonZero : lastNonZero);
            illumProfNew = IlluminanceProfile(Interp2D(x, illumReqNew, 'cubic'));
            intenProfNew = IntensityProfile(intenProf.Func, alphaMin, alphaMax);
            
            % compute norm const
            normConst = AnalyticalEngine.computeNormConstInten2Illum(intenProfNew, illumProfNew);
            
            % compute ray-corresponding function
            [psi, ro2] = ode45(@AnalyticalEngine.diffeqnReflCorrFuncUpInten2Illum, psiMir, (illumProfNew.RMax)^2, [], intenProfNew, illumProfNew, normConst);
            rOut = sqrt(abs(ro2));            
            % remove bad values
            rOut(isinf(rOut)) = NaN;
            [~, iMin] = min(abs(rOut));
            rOut = interp1( (psi(1:iMin) - psi(1))/(psi(iMin) - psi(1))*(alphaMax - alphaMin) + alphaMin,...
                (rOut(1:iMin) - rOut(1))/(rOut(iMin) - rOut(1))*(illumProfNew.RMin - illumProfNew.RMax) + illumProfNew.RMax,...
                (psi - psi(1))/(psi(end) - psi(1))*(alphaMax - alphaMin) + alphaMin);
        end % function
        
        function [psi, beta] = computeReflCorrFuncUpInten2Inten(intenProfIn, intenProfOut, alphaMin, alphaMax, numPoints)
            
            % unpack
            psiIn1 = intenProfIn.PsiMin;
            psiOut1 = intenProfOut.PsiMin;
            psiOut2 = intenProfOut.PsiMax;
            
            % check if out func is bad
            res = fminbnd(@(r)intenProfOut.eval(r), psiOut1, psiOut2);
            if intenProfOut.eval(res) < 1e-6
                ME = MException('CalculateZA:BadEvaluationFunction', ...
                    'Zero or negative irradiance function value is not appropriate');
                throw(ME);
            end
            
            % preparing reqIntenNew
            dpsi = (psiOut2 - psiOut1) / (numPoints - 1);
            psi = psiOut1 : dpsi : psiOut2;
            dpsiMir = (alphaMax - alphaMin) / (numPoints - 1);
            psiMir = alphaMin : dpsiMir : alphaMax;
            dpsiDir = (alphaMin - psiIn1) / (numPoints - 1);
            psiDir = psiIn1 : dpsiDir : alphaMin;
            intenVal = IntensityAxisym(intenProfOut).scaleEval(0, psi, 1);
            intenDir = IntensityAxisym(intenProfIn).scaleEval(0, psiDir, 1);    
            intenDirInterp = Interp2D(psiDir, intenDir, 'cubic');
            intenDir = intenDirInterp.eval(psi);
            intenDir(isnan(intenDir)) = 0;    
            intenReqNew = intenVal - intenDir;
            intenReqNew(intenReqNew < 0) = 0;
            
            intenProfOutNew = IntensityProfile(Interp2D(psi, intenReqNew, 'cubic'));
            intenProfInNew = IntensityProfile(intenProfIn.Func, alphaMin, alphaMax);
            
            % compute norm const
            normConst = AnalyticalEngine.computeNormConstInten2Inten(intenProfInNew, intenProfOutNew);
            
            % compute ray-corresponding function
            [psi, cosBeta] = ode45(@AnalyticalEngine.diffeqnReflCorrFuncUpInten2Inten, psiMir, cos(intenProfOutNew.PsiMax), [], intenProfInNew, intenProfOutNew, normConst);
            firstBadCosBeta = find(abs(cosBeta) > 1 | isinf(cosBeta) | isnan(cosBeta), 1, 'first');
            if ~isempty(firstBadCosBeta) 
                firstBadCosBeta = firstBadCosBeta - 1;
                cosBeta = interp1( (psi(1:firstBadCosBeta) - psi(1))/(psi(firstBadCosBeta) - psi(1))*(alphaMax - alphaMin) + alphaMin,...
                    (cosBeta(1:firstBadCosBeta) - cosBeta(1))/(cosBeta(firstBadCosBeta) - cosBeta(1))*(cos(intenProfOutNew.PsiMin) - cos(intenProfOutNew.PsiMax)) + cos(intenProfOutNew.PsiMax),...
                    (psi - psi(1))/(psi(end) - psi(1))*(alphaMax - alphaMin) + alphaMin );
            end
            beta = acos(cosBeta);
            
        end % function
               
    end % methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PRIVATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Static = true, Access = private)
        
        function [roErr, psi, r] = getProfileInten2Illum(eff, intenProf, illumProf, R0, nRef, z0, ...
                                                                                normConst, numPoints)
            % unpack parameters
            psiMin = intenProf.PsiMin;
            psiMax = intenProf.PsiMax;
            Rad1 = illumProf.RMin;
            Rad2 = illumProf.RMax;
            
            % get lens profile
            dpsi = abs(psiMax - psiMin) / (numPoints - 1);
            psi = psiMin : dpsi : psiMax;
            [psi, t] = ode45(@AnalyticalEngine.diffeqnProfInten2Illum, psi, [R0 Rad1^2], [], intenProf, illumProf, nRef, z0, normConst, eff);
            r = t(:,1);
            ro = sqrt(2 * t(:,2));
            
            % get ro error
            if ~isreal(ro(end)) || isnan(ro(end))
                ro(end) = 1e15;
            end % if
            roErr = ro(end) - Rad2;
        end % function

        function [bettaErr, psi, r] = getProfileInten2Inten(eff, intenProfIn, intenProfOut, R0, nRef, normConst, numPoints)
            % unpack parameters
            psiMin = intenProfIn.PsiMin;
            psiMax = intenProfIn.PsiMax;
            bettaMin = intenProfOut.PsiMin;
            bettaMax = intenProfOut.PsiMax;

            % get lens profile
            dpsi = abs(psiMax - psiMin) / (numPoints - 1);
            psi = psiMin : dpsi : psiMax;
            try
                [psi, lp] = ode45(@AnalyticalEngine.diffeqnProfInten2Inten, psi, [cos(bettaMin) R0], [], intenProfIn, intenProfOut, nRef, normConst, eff);
            catch
                bettaErr = 2*pi;
                psi = [];
                r = [];
                return;
            end
            r = lp(:,2);
            
            % get betta err
            cosBetta = lp(:,1);
            betta = acos(cosBetta);
            if ~isreal(betta(end)) || isnan(betta(end))
                betta(end) = 1e15;
            end % if
            bettaErr = betta(end) - bettaMax;
        end % function

        function drRo = diffeqnProfInten2Illum(psi, rRo, intenProf, illumProf, nRef, z0, normConst, eff)
            % function beta(psi)
            beta = atan((rRo(1) * sin(psi) - sqrt(2 * rRo(2)) ) / ...
                (z0 - rRo(1) * cos(psi)) );
            % fresnel losses
            a0 = [sin(psi) cos(psi)];
            a1 = [-sin(beta) cos(beta)];
            vn = -(a1 - nRef * a0) / sqrt(1 + nRef^2 - 2 * nRef * dot(a1,a0));
            ang_a0_vn = acos(dot(a0, vn));
            T_coeff = fresnelT(ang_a0_vn, nRef, 1);
            
            % THE FIRST DIFFERENTIAL EQUATION
            drRo(1) = -rRo(1) * sin(psi + beta) / (nRef - cos(psi + beta));
            
            % THE SECOND DIFFERENTIAL EQUATION
            % get E
            r = sqrt(rRo(2) * 2);
            E = eff * normConst * illumProf.eval(r)/1e6;
            % get I0
            I0 = intenProf.eval(psi);
            % the second differential equation
            drRo(2) = I0 * sin(psi) * T_coeff ./ E;
            
            drRo = drRo';
        end % function
        
        function dcosBettaR = diffeqnProfInten2Inten(psi, cosBettaR, intenProfIn, intenProfOut, nRef, normConst, eff)
            cosBetta = cosBettaR(1);
            betta = acos(cosBetta);
            r = cosBettaR(2);
            IIn = intenProfIn.eval(psi);
            IOut = intenProfOut.eval(betta);
            
            % set of 2 differential equations
            % fresnel losses
            a0 = [sin(psi) cos(psi)];
            a1 = [sin(betta) cos(betta)];
            vn = -(a1 - nRef * a0) / sqrt(1 + nRef^2 - 2 * nRef * dot(a1,a0));
            ang_a0_vn = acos(dot(a0, vn));
            T_coeff = fresnelT(ang_a0_vn, nRef, 1);
            
            % function beta(psi)
            dcosBettaR(1) = - IIn * sin(psi) * T_coeff / normConst / IOut / eff ;
            % THE FIRST DIFFERENTIAL EQUATION
            dcosBettaR(2) = r * sin(betta - psi) / (nRef - cos(betta - psi));
            dcosBettaR = dcosBettaR';
        end % function
        
        function dzdr = diffeqnProfColIllum2Inten(r, z, rIn, psiOut, nRef)
            psiOutCur = interp1(rIn, psiOut, r);
            dzdr = sin(psiOutCur) ./ (nRef - cos(psiOutCur));
        end
        
        function dzdr = diffeqnProfColIllum2Illum(r, z, r0, r_out, z0, nRef)
            r_out_cur = interp1(r0, r_out, r);
            a1x = r_out_cur - r;
            a1z = z0 - z;
            a1 = sqrt(a1x .^ 2 + a1z .^ 2);
            a1x = a1x ./ a1;
            a1z = a1z ./ a1;
            dzdr = a1x ./ (nRef - a1z);
        end
        
        function dlGamR = diffeqn2ProfInten2Illum(phi, lGamR, n1, n2, f, coeffs, teta)
            % Ќќ¬џ…  ќƒ
            % teta = 0 - работает только внутренн€€ поверхность ... 1 - только внешн€€, 1 пока что не работает!!!
            l = lGamR(1);
            gam = lGamR(2);
            R = lGamR(3);
            beta = (gam - teta * phi) / (1 - teta);            
            % вычисление dR_dphi
            dR = -R * sin(phi - gam) / (n1/n2 - cos(phi - gam));            
            % вычисление ro, dro
            ro = polyval(coeffs, phi);
            dro = polyval(polyder(coeffs), phi);            
            % вычисление u и v
            u = ro - (R * sin(phi) + l * sin(gam));
            v = f - (R * cos(phi) + l * cos(gam));            
            % вычисление h1, h2, h3
            h1 = (v * sin(gam) - u * cos(gam)) / v^2;
            h2 = 1 / (1 - teta) / (cos(beta))^2 + l * (v * cos(gam) + u * sin(gam)) / v^2;
            h3 = teta / (1 - teta) / (cos(beta))^2 + ...
                (v * dro - R * (v * cos(phi) + u * sin(phi)) + ...
                dR * (u * cos(phi) - v * sin(phi))) / v^2;            
            % вычисление g2, g3
            g2 = l * sin(gam - beta) / (n2 - cos(phi - gam));
            g3 = (dR * (cos(phi - beta) - n1) - R * sin(phi - beta)) / ...
                (n2 - cos(phi - gam));            
            % вычисление правых частей диффуров
            dgam = (h3 - h1 * g3) / (h2 - h1 * g2);
            dl = (h3 * g2 - h2 * g3) / (h1 * g2 - h2);            
            %%%%%%% Ќ≈ ”ƒјЋя“№ - “”“ ¬—≈ ѕќЌя“Ќќ!
            % % % l = l_gam_R(1);
            % % % gam = l_gam_R(2);
            % % % R = l_gam_R(3);
            % % %
            % % % % вычисление dR_dphi
            % % % dR = -R * sin(phi - gam) / (n1/n2 - cos(phi - gam));
            % % %
            % % % % вычисление ro, dro
            % % % ro = polyval(coeffs, phi);
            % % % dro = polyval(polyder(coeffs), phi);
            % % %
            % % % % вычисление u и v
            % % % u = ro - (R * sin(phi) + l * sin(gam));
            % % % v = f - (R * cos(phi) + l * cos(gam));
            % % %
            % % % % вычисление h1, h2, h3
            % % % h1 = (v * sin(gam) - u * cos(gam)) / v^2;
            % % % h2 = 2 / (cos(2 * gam - phi))^2 + l * (v * cos(gam) + u * sin(gam)) / v^2;
            % % % h3 = 1 / (cos(2 * gam - phi))^2 + (v * dro - R * (v * cos(phi) + u * sin(phi)) + ...
            % % %     dR * (u * cos(phi) - v * sin(phi))) / v^2;
            % % %
            % % % % вычисление g2, g3
            % % % g2 = l * sin(phi - gam) / (n2 - cos(phi - gam));
            % % % g3 = (dR * (cos(2*(phi - gam)) - n1) - R * sin(2*(phi - gam))) / ...
            % % %     (n2 - cos(phi - gam));
            % % %
            % % % % вычисление правых частей диффуров
            % % % dgam = (h3 - h1 * g3) / (h2 - h1 * g2);
            % % % dl = (h3 * g2 - h2 * g3) / (h1 * g2 - h2);            
            dlGamR = [dl dgam dR]';
        end
        
        function drInndpsi = diffeqn2ProfInten2Inten1(phi, R, n1, n2, coeffsgam)
            gam = polyval(coeffsgam, phi);
            drInndpsi = -R * sin(phi - gam) / (n1/n2 - cos(phi - gam));
        end
        
        function dlOutdpsi = diffeqn2ProfInten2Inten2(phi, l,  n1, n2, coeffsbetta, coeffsgam, coeffsR, coeffsRDer, coeffsgamDer)
            R = polyval(coeffsR, phi);
            betta = polyval(coeffsbetta, phi);
            gam = polyval(coeffsgam, phi);
            dR_dphi = polyval(coeffsRDer, phi);
            dgam_dphi = polyval(coeffsgamDer, phi);
            dlOutdpsi = (dR_dphi .* cos(phi - betta) - R .* sin(phi - betta) - dgam_dphi .* ...
                l .* sin(gam - betta) - n1 * dR_dphi) ./ ...
                (n2 - cos(phi - gam));
        end
        
        function drdpsi = diffeqnReflProfUpInten2Illum(psi, r, z0, coeffs)
            ro = polyval(coeffs, psi);
            beta = atan((r * sin(psi) - ro) / (z0 - r * cos(psi)));
            drdpsi = - r .* tan((pi - psi - beta)/2);
        end % function
        
        function drdpsi = diffeqnReflProfDownInten2Illum(psi, r, z0, coeffs)
            ro = polyval(coeffs, psi);
            beta = atan((ro - r .* sin(psi)) ./ (z0 + r .* cos(psi)));
            drdpsi = r .* tan((psi + beta)/2);
        end % function
        
        function drdpsi = diffeqnReflProfUpInten2Inten(psi, r, coeffs)
            beta = polyval(coeffs, psi);
            drdpsi = - r .* tan((pi - psi - beta)/2);
        end % function
        
        function drdpsi = diffeqnReflProfDownInten2Inten(psi, r, coeffs)
            beta = polyval(coeffs, psi);
            drdpsi = r .* tan((psi + beta)/2);
        end % function
        
        function drAdpsi = diffeqnTIRProfInten2Inten1(betaCur, r, beta0, betaA, nRef)
            betaACur = interp1(beta0, betaA, betaCur);
            drAdpsi = r * sin(betaACur - betaCur) / (1/nRef - cos(betaACur - betaCur));
        end
        
        function drCdpsi = diffeqnTIRProfInten2Inten2(phi, l, coeffsRB, coeffsGamma, coeffsRBdPhi,...
                                                        coeffsGammadPhi, coeffsS1x, coeffsS1z, nRef)
            rB = polyval(coeffsRB, phi);
            gamma = polyval(coeffsGamma, phi);
            rBdPhi = polyval(coeffsRBdPhi, phi);
            gammadPhi = polyval(coeffsGammadPhi, phi);
            s1x = polyval(coeffsS1x, phi) * nRef;
            s1z = polyval(coeffsS1z, phi) * nRef;
            
            drCdpsi = (s1x * (rBdPhi * cos(phi) - rB * sin(phi) - l * sin(gamma) *gammadPhi) + ...
                s1z * (rBdPhi * sin(phi) + rB * cos(phi) + l * cos(gamma)*gammadPhi) - rBdPhi) / ...
                (nRef - s1x * cos(gamma) - s1z * sin(gamma));
        end
        
        function drAdphi = diffeqnTIRProfInten2Illum1(betaCur, r, beta0, betaA, nRef)
            betaACur = interp1(beta0, betaA, betaCur);
            drAdphi = r * sin(betaACur - betaCur) / (1/nRef - cos(betaACur - betaCur));
        end
        
        function drCdphi = diffeqnTIRProfInten2Illum2(phi, l, coeffsRB, coeffsGamma, coeffsRBdPhi,...
                                                coeffsGammadPhi, coeffsS1x, coeffsS1z, nRef)
            rB = polyval(coeffsRB, phi);
            gamma = polyval(coeffsGamma, phi);
            rBdPhi = polyval(coeffsRBdPhi, phi);
            gammadPhi = polyval(coeffsGammadPhi, phi);
            s1x = polyval(coeffsS1x, phi) * nRef;
            s1z = polyval(coeffsS1z, phi) * nRef;
            
            drCdphi = (s1x * (rBdPhi * cos(phi) - rB * sin(phi) - l * sin(gamma) *gammadPhi) + ...
                s1z * (rBdPhi * sin(phi) + rB * cos(phi) + l * cos(gamma)*gammadPhi) - rBdPhi) / ...
                (nRef - s1x * cos(gamma) - s1z * sin(gamma));
        end
        
        function dCosPsidR = diffeqnCorrFuncIllum2Inten(r, cosPsi, illumProf, intenProf, normConst)
            psi = acos(cosPsi);
            Ival = intenProf.eval(psi);
            illumVal = illumProf.eval(r);
            dCosPsidR = - illumVal .* r ./ Ival / normConst /1e6;
        end
        
        function drOut2dr = diffeqnCorrFuncIllum2Illum(r, r_out, illumProfIn, illumProfOut, normConst)
            illumOutVal = illumProfOut.eval(sqrt(abs(r_out)));
            illumInVal = illumProfIn.eval(r);
            drOut2dr = 2 * illumInVal .* r ./ illumOutVal / normConst;
        end
        
        function drOutdPsi = diffeqnCorrFuncInten2Illum(psi, ro2, intenProf, illumProf, normConst)
            ro = sqrt(abs(ro2));
            irradVal = illumProf.eval(ro);
            intenVal = intenProf.eval(psi);
            drOutdPsi = 2 * intenVal .* sin(psi) / normConst ./ irradVal * 1e6;
        end
        
        function dpsidcosBetta = diffeqnCorrFuncInten2Inten(psi, cosBetta, intenProfIn, intenProfOut, normConst)
            betta = acos(cosBetta);
            IBeta = intenProfOut.eval(betta);
            intenInVal = intenProfIn.eval(psi);
            dpsidcosBetta = - intenInVal * sin(psi) / normConst / IBeta;
        end
        
        function prof = getColProfileCart(r0, psiBorders, nRef, numPoints)
            if length(psiBorders) == 1
                psiMin = 0;
                psiMax = psiBorders;
            else
                psiMin = psiBorders(1);
                psiMax = psiBorders(2);
            end % if
            dpsi = (psiMax - psiMin) / (numPoints - 1);
            psi = (psiMin : dpsi : psiMax)';
            r = r0 * (1 - nRef * cos(psiMin)) ./ (1 - nRef * cos(psi));
            x = r .* sin(psi); z = r .* cos(psi);
            prof = [x z];
        end
        
        function drodpsi = diffeqnReflCorrFuncUpInten2Illum(psi, ro2, intenProf, illumProf, normConst)
            ro = sqrt(abs(ro2));
            illumVal = illumProf.eval(ro);
            intenVal = intenProf.eval(psi);
            drodpsi = -2 * intenVal .* sin(psi) / normConst ./ illumVal * 1e6;            
        end % function
        
        function dpsidcosBetta = diffeqnReflCorrFuncUpInten2Inten(psi, cosBetta, intenProfIn, intenProfOut, normConst)
            betta = abs(acos(cosBetta));
            IBeta = intenProfOut.eval(betta);
            intenInVal = intenProfIn.eval(psi);
            dpsidcosBetta = intenInVal * sin(psi) / normConst / IBeta;
        end % function
        
    end % methods
    
end % classdef