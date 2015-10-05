function [p, ixSt]  = plotTimeSeries(stTS, xlim)
    if (~exist('xlim','var'))
        bx = 0;
        xlim = [];
    else
        bx = 1;
    end
    yPlots = 4;
    xPlots = 1;
    ixPlot = 0;
    bMovingAverage = 0;
    bSegmentedIntensity = 1;
    bNumberFlux = 1;

    sMarker = 'o';
    nMarkerSize = 5;
    nLineWidth = 2;

    %%%
    % Plot the oxygen profile.
    if 1
        ixPlot = ixPlot + 1;
        subplot(yPlots, xPlots, ixPlot);
        p = plotOxygenProfile(stTS);
        %set(p, 'marker', sMarker, 'markersize', nMarkerSize, 'linewidth', nLineWidth);
        if (~isempty(stTS(1).intensity))
            hold on;
%            p = plotMedianIntensity(stTS, [min([stTS(:).oxygen]) max([stTS(:).oxygen])], xlim, bMovingAverage, bSegmentedIntensity);

            p = plotMedianIntensity(stTS, prctile([stTS(:).oxygen], [0 100]), xlim, bMovingAverage, bSegmentedIntensity);

            set(p, 'marker', sMarker, 'markersize', nMarkerSize, 'linewidth', nLineWidth);

            hold off;
            %set(p, 'color', 'r');
        end
        if (bx)
            set(gca, 'xlim', xlim);
        end
    end

    if 1
        %%%
        % Plot the global pressure drop from the inlet to the outlet (atmospheric).
        ixPlot = ixPlot + 1;
        subplot(yPlots, xPlots, ixPlot);
        p = plotAppliedPressure(stTS, xlim);
        if (bx)
            set(gca, 'xlim', xlim);
        end
    end

    %%%
    % Plot the median velocity.
    ixPlot = ixPlot + 1;
    subplot(yPlots, xPlots, ixPlot);
    [p, ixSt] = plotMedianVelocity(stTS, xlim, bMovingAverage);
    if (bx)
        set(gca, 'xlim', xlim);
    end

    %%%
    % Plot the effective viscosity.
    ixPlot = ixPlot + 1;
    subplot(yPlots, xPlots, ixPlot);
    [p, ixSt] = plotApparentViscosity(stTS, xlim, bMovingAverage);
    if (bx)
        set(gca, 'xlim', xlim);
    end

    if 0
        ixPlot = ixPlot + 1;
        subplot(yPlots, xPlots, ixPlot);
        p = plotMedianHematocrit(stTS);
        set(p, 'marker', sMarker, 'markersize', nMarkerSize, 'linewidth', 2);
        if (bx)
            set(gca, 'xlim', xlim);
        end
    end

    if 0
        ixPlot = ixPlot + 1;
        subplot(yPlots, xPlots, ixPlot);
        p = plotNetFlow(stTS, bNumberFlux);
        if (bx)
            set(gca, 'xlim', xlim);
        end
        ixPlot = ixPlot + 1;
        subplot(yPlots, xPlots, ixPlot);
        p = plotMedianVelocitySpread(stTS);
        if (bx)
            set(gca, 'xlim', xlim);
        end

        ixPlot = ixPlot + 1;
        subplot(yPlots, xPlots, ixPlot);
        p = plotHematocritSpread(stTS);
        if (bx)
            set(gca, 'xlim', xlim);
        end
    end

end

function p = plotNetFlow(stTS, bNumberFlux)
    vTime = [stTS(:).time]/60;
    vNetFlow = NaN*ones(numel(stTS), 1);
    for ix = 1:numel(stTS)
        if (bNumberFlux)
            vNetFlow(ix) = abs(stTS(ix).velocity(1)*stTS(ix).cellNumber(1) - stTS(ix).velocity(end - 2)*stTS(ix).cellNumber(end - 2));
        else
            vNetFlow(ix) = abs(stTS(ix).velocity(1)*stTS(ix).hematocrit(1) - stTS(ix).velocity(end - 2)*stTS(ix).hematocrit(end - 2));
        end
    end
    [vBinTime, vBinNetFlow] = binMeasurements(vTime, vNetFlow, 120);
    %p = plot(vBinTime, vBinNetFlow, 'linestyle', 'none', 'marker', '.', 'markersize', 8);
    p = plot(vTime, vNetFlow, 'linestyle', 'none', 'marker', '.', 'markersize', 8);
    xlabel('minutes');
    if (bNumberFlux)
        ylabel('net flow (#*m/s)');
    else
        ylabel('net flow (%*m/s)');
    end
end

% Sample the data.
function [ox, oy] = binMeasurements(ix, iy, n)
    n = min(n, length(ix));
    fStep = length(iy)/n;
    if (fStep == 1)
        ox = ix;
        oy = iy;
    else
        ox = NaN*ones(n,1);
        oy = ox;
        for jx = 1:n
            kx2 = floor(fStep*jx);
            kx1 = max(1, ceil(fStep*(jx - 1)));
            ox(jx) = mean(ix(kx1:kx2));
            ty = iy(kx1:kx2);
            oy(jx) = mean(ty(isfinite(ty)));
        end
    end
end

function p = plotHematocritSpread(stTS)
    vTime = [stTS(:).time]/60;
    vHematocritSpread = NaN*ones(numel(stTS), 1);
    for ix = 1:numel(stTS)
        vHematocritSpread(ix) = stTS(ix).hematocrit(end);
    end
    p = plot(vTime, vHematocritSpread, 'linestyle', 'none', 'marker', '.', 'markersize', 8);
    xlabel('minutes');
    ylabel('hematocrit spread (%)');
end

function p = plotMedianHematocrit(stTS)
    vTime = [stTS(:).time]/60;
    vMedianHematocrit = NaN*ones(numel(stTS), 1);
    for ix = 1:numel(stTS)
        vMedianHematocrit(ix) = stTS(ix).hematocrit(end - 1);
    end
    p = plot(vTime, vMedianHematocrit, 'linestyle', 'none', 'marker', '.', 'markersize', 8);
    xlabel('minutes');
    ylabel('median hematocrit (%)');
end

function p = plotMedianIntensity(stTS, vScale, xlim, bMovingAverage, bSegmentedIntensity)
    bMovingAverage = 1;
    vTime = [stTS(:).time]/60;
    vMedianIntensity = NaN*ones(numel(stTS), 1);
    if (bSegmentedIntensity)
        ixCol = 2;
    else
        ixCol = 1;
    end

    if (isfield(stTS(1), 'index'))
        nMaxIndex = max([stTS(:).index]);
        vMedianIntensity = NaN*ones(nMaxIndex, length(stTS));
        for ix = 1:numel(stTS)
            vMedianIntensity(stTS(ix).index, ix) = stTS(ix).intensity(end - 1, ixCol);
        end
    else
        for ix = 1:numel(stTS)
            vMedianIntensity(ix) = stTS(ix).intensity(end - 1, ixCol);
        end
    end

    if (~isempty(xlim))
        ix = intersect(find(vTime > xlim(1)), find(vTime < xlim(2)));
    else
        ix = 1:max(size(vMedianIntensity));
    end
    nMax = max(vMedianIntensity(1, ix));
    nMin = min(vMedianIntensity(1, ix));

    bScale = 0
    if (bScale)
        if (exist('vScale','var'))
            vMedianIntensity = vMedianIntensity - nMin;
            %vMedianIntensity = vMedianIntensity - 0;
            vMedianIntensity = vMedianIntensity/(nMax - nMin);
            %vMedianIntensity = vMedianIntensity/max(vMedianIntensity(isfinite(vMedianIntensity)));
            %vMedianIntensity = vMedianIntensity * (vScale(2) - vScale(1) + 10);
            %vMedianIntensity = vMedianIntensity + vScale(1) - 5;
            vMedianIntensity = vMedianIntensity * (vScale(2) - vScale(1));
            vMedianIntensity = vMedianIntensity + vScale(1);
        end
    end
    
    if (bMovingAverage)
        for ix = 1:size(vMedianIntensity, 2)
            vMedMedIntensity(ix) = median(vMedianIntensity(isfinite(vMedianIntensity(:, ix)), ix));
        end
        n = 1;
        w = conv([1:n]/sum(1:n), vMedMedIntensity);
        w = smooth(vMedMedIntensity, n);
        p = plot(vTime, w(1:length(vTime)), 'linestyle', 'none', 'marker', '.', 'markersize', 8);
    else
        p = plot(vTime, vMedianIntensity, 'linestyle', 'none', 'marker', '.', 'markersize', 8);
    end
    xlabel('minutes');
    ylabel('median pixel intensity');
end


function [p, ixSt] = plotApparentViscosity(stTS, xlim, bMovingAverage)
    p = [];
    vTime = [stTS(:).time]/60;

    % Crop to time limits.
    if (~isempty(xlim))
        ixSt = intersect(find(vTime >= xlim(1)), find(vTime <= xlim(2)));
        stTS = stTS(ixSt);
        vTime = [stTS(:).time]/60;
    else
        ixSt = [1:size(stTS,1)];
    end
    vMedianVelocity = NaN*ones(numel(stTS), 1);
    
    % Initialize cell arrays for velocities.
    nv = size(stTS(1).velocity, 1);
    if (isfield(stTS(1), 'index'))
        nMaxIndex = max([stTS(:).index]);
        csMedianVelocity = cell(nMaxIndex, 1);
        csFilteredVelocity = cell(nMaxIndex, 1);
        csApparentViscosity = cell(nMaxIndex, 1);
        for ix = 1:nMaxIndex
            csMedianVelocity{ix} = vMedianVelocity;
            csFilteredVelocity{ix} = vMedianVelocity;
            csApparentViscosity{ix} = vMedianVelocity;
        end
        nv = 1;
    else
        csMedianVelocity = cell(nv, 1);
        csFilteredVelocity = cell(nv, 1);
    end
    
    % Collect all the median velocities in a cell array indexed by (time,
    % index).
    for ix = 1:numel(stTS)
        %vMedianVelocity(ix) = stTS(ix).velocity(end - 1);
        for jx = 1:nv
            % The last velocity row is the median.
            if (jx <= size(stTS(ix).velocity, 1))
                if (size(stTS(ix).velocity, 2) > 1)
                    nMedianVelocity = stTS(ix).velocity(jx, 1);
                    nVelocityIQR = stTS(ix).velocity(jx, 2);
                else
                    nMedianVelocity = stTS(ix).velocity(end - 1);
                    nVelocityIQR = stTS(ix).velocity(end);
                    if (isnan(nMedianVelocity))
                        ixFinite = find(isfinite(stTS(ix).velocity([1:(end - 2)])));
                    end
                end

                % Ignore if the velocity is greater than 50 pixels/s and the
                % variation is greater than the median velocity multiplied by the noiselevel.
                kNoiseLevel = 1;
                nVelocityThreshold = 200;
                if (abs(nMedianVelocity) > nVelocityThreshold)
                    if (kNoiseLevel*abs(nMedianVelocity) < nVelocityIQR)
                        nFilteredVelocity = nMedianVelocity;
                        nMedianVelocity = NaN;
                    else
                        nFilteredVelocity = NaN;
                    end
                else
                    if (nVelocityIQR > nVelocityThreshold*kNoiseLevel)
                        nMedianVelocity = NaN;
                        nFilteredVelocity = nMedianVelocity;
                    else
                        nFilteredVelocity = NaN;
                    end
                end
                if (isfield(stTS(1), 'index'))
                    kx = stTS(ix).index;
                else
                    kx = jx;
                end
                csMedianVelocity{kx}(ix) = nMedianVelocity;
                csFilteredVelocity{kx}(ix) = nFilteredVelocity;
                csApparentViscosity{kx}(ix) = nMedianVelocity/stTS(ix).pressure;
            else
                csMedianVelocity{jx}(ix) = NaN;
                csFilteredVelocity{jx}(ix) = NaN;
                csApparentViscosity{kx}(ix) = NaN;
            end
        end
        vMedianVelocity(ix) = stTS(ix).velocity(1, 1);
    end

    % Plot only a moving average?
    if (bMovingAverage)
        vWeights = 1.5.^(1:8);
        w = conv(vWeights/sum(vWeights), vMedianVelocity);
        p = plot(vTime, w(1:length(vTime)), 'linestyle', 'none', 'marker', '.', 'markersize', 8);
    else
        hold on
        if (isfield(stTS(1), 'index'))
            nLines = nMaxIndex;
        else
            nLines = nv -1;
        end
        cmap = lines(nLines);
        for jx = 1:nLines
            [vBinTime, vBinVelocity] = binMeasurements(vTime, csApparentViscosity{jx}, round(length(stTS)/1));
            p = plot(vBinTime, abs(vBinVelocity(:)), 'color', cmap(jx, :), 'marker', 'o', 'markersize', 6, 'linewidth', 1.5, 'linestyle', 'none');
            [vBinTime, vBinVelocity] = binMeasurements(vTime, csApparentViscosity{jx}, round(length(stTS)/1));
            p = plot(vBinTime, abs(vBinVelocity(:)), 'color', cmap(jx, :), 'marker', 'x', 'markersize', 7, 'linewidth', 1.5, 'linestyle', 'none');
        end
        pl = legend('show','location','west');
        set(pl, 'fontsize', 6);
    end
    xlabel('minutes');
    ylabel('1/apparent viscosity');
end

function [p, ixSt] = plotMedianVelocity(stTS, xlim, bMovingAverage)
    p = [];
    vTime = [stTS(:).time]/60;

    % Crop to time limits.
    if (~isempty(xlim))
        ixSt = intersect(find(vTime >= xlim(1)), find(vTime <= xlim(2)));
        stTS = stTS(ixSt);
        vTime = [stTS(:).time]/60;
    else
        ixSt = [1:size(stTS,1)];
    end
    vMedianVelocity = NaN*ones(numel(stTS), 1);
    
    % Initialize cell arrays for velocities.
    nv = size(stTS(1).velocity, 1);
    if (isfield(stTS(1), 'index'))
        nMaxIndex = max([stTS(:).index]);
        csMedianVelocity = cell(nMaxIndex, 1);
        csFilteredVelocity = cell(nMaxIndex, 1);
        csApparentViscosity = cell(nMaxIndex, 1);
        for ix = 1:nMaxIndex
            csMedianVelocity{ix} = vMedianVelocity;
            csFilteredVelocity{ix} = vMedianVelocity;
            csApparentViscosity{ix} = vMedianVelocity;
        end
        nv = 1;
    else
        csMedianVelocity = cell(nv, 1);
        csFilteredVelocity = cell(nv, 1);
    end
    
    % Collect all the median velocities in a cell array indexed by (time,
    % index).
    for ix = 1:numel(stTS)
        %vMedianVelocity(ix) = stTS(ix).velocity(end - 1);
        for jx = 1:nv
            % The last velocity row is the median.
            if (jx <= size(stTS(ix).velocity, 1))
                if (size(stTS(ix).velocity, 2) > 1)
                    nMedianVelocity = stTS(ix).velocity(jx, 1);
                    nVelocityIQR = stTS(ix).velocity(jx, 2);
                else
                    nMedianVelocity = stTS(ix).velocity(end - 1);
                    nVelocityIQR = stTS(ix).velocity(end);
                    if (isnan(nMedianVelocity))
                        ixFinite = find(isfinite(stTS(ix).velocity([1:(end - 2)])));
                    end
                end

                % Ignore if the velocity is greater than 50 pixels/s and the
                % variation is greater than the median velocity multiplied by the noiselevel.
                kNoiseLevel = 1;
                nVelocityThreshold = 200;
                if (abs(nMedianVelocity) > nVelocityThreshold)
                    if (kNoiseLevel*abs(nMedianVelocity) < nVelocityIQR)
                        nFilteredVelocity = nMedianVelocity;
                        nMedianVelocity = NaN;
                    else
                        nFilteredVelocity = NaN;
                    end
                else
                    if (nVelocityIQR > nVelocityThreshold*kNoiseLevel)
                        nMedianVelocity = NaN;
                        nFilteredVelocity = nMedianVelocity;
                    else
                        nFilteredVelocity = NaN;
                    end
                end
                if (isfield(stTS(1), 'index'))
                    kx = stTS(ix).index;
                else
                    kx = jx;
                end
                csMedianVelocity{kx}(ix) = nMedianVelocity;
                csFilteredVelocity{kx}(ix) = nFilteredVelocity;
                csApparentViscosity{kx}(ix) = nMedianVelocity/stTS(ix).pressure;
            else
                csMedianVelocity{jx}(ix) = NaN;
                csFilteredVelocity{jx}(ix) = NaN;
                csApparentViscosity{kx}(ix) = NaN;
            end
        end
        vMedianVelocity(ix) = stTS(ix).velocity(1, 1);
    end

    % Plot only a moving average?
    if (bMovingAverage)
        vWeights = 1.5.^(1:8);
        w = conv(vWeights/sum(vWeights), vMedianVelocity);
        p = plot(vTime, w(1:length(vTime)), 'linestyle', 'none', 'marker', '.', 'markersize', 8);
    else
        hold on
        if (isfield(stTS(1), 'index'))
            nLines = nMaxIndex;
        else
            nLines = nv -1;
        end
        cmap = lines(nLines);
        for jx = 1:nLines
            [vBinTime, vBinVelocity] = binMeasurements(vTime, csMedianVelocity{jx}, round(length(stTS)/1));
            %p = plot(vBinTime, abs(vBinVelocity(:)), 'color', cmap(jx, :), 'marker', 'o', 'markersize', 6, 'linewidth', 1.5, 'linestyle', 'none');
            p = plot(vBinTime, (vBinVelocity(:)), 'color', cmap(jx, :), 'marker', 'o', 'markersize', 6, 'linewidth', 1.5, 'linestyle', 'none');
            [vBinTime, vBinVelocity] = binMeasurements(vTime, csFilteredVelocity{jx}, round(length(stTS)/1));
            p = plot(vBinTime, abs(vBinVelocity(:)), 'color', cmap(jx, :), 'marker', 'x', 'markersize', 7, 'linewidth', 1.5, 'linestyle', 'none');
        end
        pl = legend('show','location','west');
        set(pl, 'fontsize', 6);
    end
    xlabel('minutes');
    ylabel('velocity (pixels/s)');
end

function p = plotMedianVelocitySpread(stTS)
    vTime = [stTS(:).time]/60;
    vVelocitySpread = NaN*ones(numel(stTS), 1);
    for ix = 1:numel(stTS)
        vVelocitySpread(ix) = stTS(ix).velocity(end)/stTS(ix).velocity(end - 1);
    end
    p = plot(vTime, vVelocitySpread, 'linestyle', 'none', 'marker', '.', 'markersize', 8);
    set(gca,'ylim',[0 1]);
    xlabel('minutes');
    ylabel('velocity spread (%)');
end

function p = plotAppliedPressure(stTS, xlim)
    vTime = [stTS(:).time]/60;
    vPressure = [stTS(:).pressure];

    if (isfield(stTS(1), 'index'))
        sLinestyle = 'none';
    else
        sLinestyle = '-';
    end
    p = plot(vTime, vPressure, 'linestyle', sLinestyle , 'marker', '.', 'markersize', 5);
    xlabel('minutes');
    ylabel('pressure (psi)');
end

function p = plotOxygenProfile(stTS)
    vTime = [stTS(:).time]/60;
    vOxygen = [stTS(:).oxygen];
    
    if (isfield(stTS(1), 'index'))
        sLinestyle = 'none';
    else
        sLinestyle = '-';
    end
    p = plot(vTime, vOxygen, 'linestyle', sLinestyle, 'marker', 'x', 'markersize', 5, 'color', 'k');

    xlabel('minutes');
    ylabel('oxygen %');
end

% Plot an array of panels to show what is happening during an experiment.
function garbage()
    nPlots = 8;
    ixPlot = 1;
    nLinewidth = 3;
    bInterp = 0;
    colInterp = 'k';
    nTS = numel(stTS);

    % Setup the time vector.
    vTime = cat(1, stTS(:).time);
    [vTime, ixTime] = sort(vTime);
    nTimeZero = min(vTime);
    kSecondsPerDay = 60*60*24;
    vTime = (vTime - nTimeZero)*kSecondsPerDay;
    vTimeInterp = union(vTime, linspace(min(vTime), max(vTime), 100));

    % Setup the oxygen vector.
    vOxygen = cat(1, stTS(:).oxygen);
    vOxygen = vOxygen(ixTime);

    % Plot the density.
    subplot(nPlots, 1, ixPlot);
    vDensity = cat(1, stTS(:).density);
    vDensity = vDensity(ixTime);
    vDensityInterp = interp1(vTime, vDensity, vTimeInterp, 'pchip');
    plot(vTime, vDensity, 'linewidth', nLinewidth);
    hold on;
    addOxygenColor(vTime, vOxygen);
    if (bInterp), plot(vTimeInterp, vDensityInterp, colInterp); end;
    hold off;
    ylabel('Density');
    ixPlot = ixPlot + 1;

    % Plot the heterogeneity.
    subplot(nPlots, 1, ixPlot);
    vHeterogeneity = cat(1, stTS(:).heterogeneity);
    vHeterogeneity = vHeterogeneity(ixTime);
    vHeterogeneityInterp = interp1(vTime, vHeterogeneity, vTimeInterp, 'pchip');
    plot(vTime, vHeterogeneity, 'linewidth', nLinewidth);
    hold on;
    addOxygenColor(vTime, vOxygen);
    if (bInterp), plot(vTimeInterp, vHeterogeneityInterp, colInterp); end;
    ylabel('N''hood Radius');
    hold off;
    ixPlot = ixPlot + 1;

    % Plot the granularity.
    subplot(nPlots, 1, ixPlot);
    vGranulometry = cat(1, stTS(:).granulometry);
    vGranulometry = vGranulometry (ixTime);
    vGranulometryInterp = interp1(vTime, vGranulometry, vTimeInterp, 'pchip');
    plot(vTime, vGranulometry, 'linewidth', nLinewidth);
    hold on;
    addOxygenColor(vTime, vOxygen);
    if (bInterp), plot(vTimeInterp, vGranulometryInterp, colInterp); end;
    ylabel('Object Radius');
    hold off;
    ixPlot = ixPlot + 1;

    % Plot the velocity.
    subplot(nPlots, 1, ixPlot);
    vVelocity = cat(1, stTS(:).velocity);
    vVelocity = vVelocity(ixTime);
    ix = ~isnan(vVelocity);
    vVelocityInterp = interp1(vTime(ix), vVelocity(ix), vTimeInterp, 'pchip');
    plot(vTime, vVelocity, 'linewidth', nLinewidth);
    hold on;
    addOxygenColor(vTime, vOxygen);
    if (bInterp), plot(vTimeInterp, vVelocityInterp, colInterp); end;
    ylabel('Velocity (um/s)');
    hold off;
    ixPlot = ixPlot + 1;

    % Plot the distance.
    subplot(nPlots, 1, ixPlot);
    vDistance = cat(1, stTS(:).distance);
    vDistance = vDistance(ixTime);
    vDistance = connectDistances(vTime, vVelocity, vDistance);
    ix = ~isnan(vDistance);
    vDistanceInterp = interp1(vTime(ix), vDistance(ix), vTimeInterp, 'pchip');
    plot(vTime, vDistance, 'linewidth', nLinewidth);
    hold on;
    addOxygenColor(vTime, vOxygen);
    if (bInterp), plot(vTimeInterp, vDistanceInterp, colInterp); end;
    ylabel('Distance (um)');
    hold off;
    ixPlot = ixPlot + 1;

    % Plot the oxygen.
    subplot(nPlots, 1, ixPlot);
    ix = ~isnan(vOxygen);
    vOxygenInterp = interp1(vTime(ix), vOxygen(ix), vTimeInterp, 'pchip');
    plot(vTime, vOxygen, 'linewidth', nLinewidth);
    hold on;
    if (bInterp), plot(vTimeInterp, vOxygenInterp, colInterp); end;
    ylabel('Oxygen (%)');
    hold off;
    ixPlot = ixPlot + 1;

    % Plot the pressure.
    subplot(nPlots, 1, ixPlot);
    vPressure = cat(1, stTS(:).pressure);
    vPressure = vPressure(ixTime);
    ix = ~isnan(vPressure);
    vPressureInterp = interp1(vTime(ix), vPressure(ix), vTimeInterp, 'pchip');
    plot(vTime, vPressure, 'linewidth', nLinewidth);
    hold on;
    addOxygenColor(vTime, vOxygen);
    if (bInterp), plot(vTimeInterp, vPressureInterp, colInterp); end;
    ylabel('Pressure (cm water)');
    hold off;
    ixPlot = ixPlot + 1;

    % Plot the video quality.
    subplot(nPlots, 1, ixPlot);
    vQuality = cat(1, stTS(:).quality);
    vVideoTime = NaN*ones(length(stTS), 1);
    for ixVideo = 1:length(stTS)
        vVideoTime(ixVideo) = stTS(ixVideo).time(1);
    end
    [vVideoTime, ixVideo] = sort(vVideoTime);
    vVideoTime = (vVideoTime - min(vVideoTime))*kSecondsPerDay;
    vQuality = vQuality(ixVideo, :);
    vQualityInterp = interp1(vVideoTime, vQuality, vTimeInterp, 'pchip');
    %    plot(vVideoTime, vQuality(:, 1), 'linewidth', nLinewidth);
    ix = ~isnan(vVelocity);
    vCV = sqrt(var(vVelocity(ix)))/mean(vVelocity(ix));
    plot(vVideoTime, vQuality(:,1).^4.*vQuality(:, 3).*vQuality(:, 4)/vCV, 'o');
    hold on;
    %    addOxygenColor(vTime, vOxygen);
    if (bInterp), plot(vTimeInterp, vQualityInterp, colInterp); end;
    ylabel('Video Quality');
    hold off;
    ixPlot = ixPlot + 1;
end

function addOxygenColor(vTime, vOxygen)
    % Add background patch objects denoting oxygen status.
    ixChange = find(diff(vOxygen)) + 1;
    ylim = get(gca, 'ylim');
    vy = [ylim(1) ylim(1) ylim(2) ylim(2)];
    xlim = get(gca, 'xlim');
    ixStart = xlim(1);
    colHi = [1 0 0];
    colLo = [0 0 1];
    nOxygenThreshold = 3;
    nAlpha = 0.25;

    for ix_ixChange = 1:(length(ixChange))
        ixEnd = vTime(ixChange(ix_ixChange));
        vx = [ixStart ixEnd ixEnd ixStart];
        %        if (vOxygen(ixChange(ix_ixChange)) > nOxygenThreshold)
        if (vOxygen(ixChange(ix_ixChange) - 1) > nOxygenThreshold)
            tCol = colHi;
        else
            tCol = colLo;
        end
        patch(vx, vy, tCol, 'facealpha', nAlpha)
        ixStart = ixEnd;
    end

    if (vOxygen(length(vOxygen)) > nOxygenThreshold)
        tCol = colHi;
    else
        tCol = colLo;
    end
    vx = [ixStart xlim(2) xlim(2) ixStart];
    patch(vx, vy, tCol, 'facealpha', nAlpha);
end

function p = plotTimeSerieslocal(mMedianTimeSeries)
    %% Plot the median time series and save PDF and fig files.
    p = figure;
    hold on
    sMarker = '.';
    % Plot the velocities.
    nPlots = 5;
    subplot(nPlots, 1, 1);
    plot(mMedianTimeSeries(:, 1), mMedianTimeSeries(:, 2),'color','r','marker',sMarker,'linestyle','-');
    % Plot a few labels showing movie indices.
    for jx = 1:2:size(mMedianTimeSeries,1)
        p = text(mMedianTimeSeries(jx,1), mMedianTimeSeries(jx,2),sprintf('%d', mMedianTimeSeries(jx,3)));
        set(p, 'fontsize', 10);
    end
    %title(dMATFiles(1).name);
    ylabel('Velocity (microns/second)');

    % Plot the frame intensity.
    subplot(nPlots, 1, 2);
    plot(mMedianTimeSeries(:, 1), mMedianTimeSeries(:, 4),'color','r','marker',sMarker,'linestyle','-');
    set(gca,'ylim',[0 255]);
    ylabel('Intensity');

    % Plot the intracellular intensity.
    subplot(nPlots, 1, 3);
    plot(mMedianTimeSeries(:, 1), mMedianTimeSeries(:, 7),'color','r','marker',sMarker,'linestyle','-');
    ylabel('Intracellular Intensity');
    xlabel('Time (s)');

    % Plot the cell density.
    subplot(nPlots, 1, 4);
    plot(mMedianTimeSeries(:, 1), mMedianTimeSeries(:, 6),'color','r','marker',sMarker,'linestyle','-');
    set(gca,'ylim',[0 1]);
    ylabel('f(\rho)');
    xlabel('Time (s)');

    % Plot the oxygen concentration.
    subplot(nPlots, 1, 5);
    plot(mMedianTimeSeries(:, 1), mMedianTimeSeries(:, 5),'color','r','marker',sMarker,'linestyle','-');
    set(gca,'ylim',[0 25]);
    ylabel('O_{2}');
    xlabel('Time (s)');

    hold off



    if (kSaveFigure)
        p = plotTimeSerieslocal(stTimeSeries);
        print(gcf,'-dpdf',['TS.' dMATFiles(1).name '.' datestr(now, 'yyyymmdd') datestr(now,'HHMMSS') '.pdf']);
        hgsave(gcf,['TS.' dMATFiles(1).name '.' datestr(now, 'yyyymmdd') datestr(now,'HHMMSS') '.fig']);
    end
end


% function vDistance = connectDistances(vTime, vVelocity, vDistance)
%     % Make the distances cumulative.
%
%     % Start at 0.
%     if (isnan(vDistance(1)))
%         vDistance(1) = 0;
%     end
%
%     ixNaN = find(isnan(vDistance));
%     for ixJump = 1:length(ixNaN)
%         nGapVelocity = mean(vVelocity(ixNaN(ixJump) + [-1 1]));
%         nGapDuration = vTime(ixNaN(ixJump)) - vTime(ixNaN(ixJump) - 1);
%         vDistance(ixNaN(ixJump)) = nGapVelocity*nGapDuration + vDistance(ixNaN(ixJump) - 1);
%         if (ixJump == length(ixNaN))
%             ixEnd = length(vDistance);
%         else
%             ixEnd = ixNaN(ixJump + 1) - 1;
%         end
%         ixStart = (ixNaN(ixJump) + 1);
%         if ~isnan(vDistance(ixNaN(ixJump)))
%             vDistance(ixStart:ixEnd) = vDistance(ixStart:ixEnd) + vDistance(ixNaN(ixJump));
%         end
%     end
%