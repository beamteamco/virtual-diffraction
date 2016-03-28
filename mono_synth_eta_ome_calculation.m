%% Forward modeling of eta-omega maps for the cubic and mono phases
% 1. Calculate eta-omega from cubic grain orientation using virtual
%    diffractometer
% 2. Superimpose synth spots on GE data. Plot it. [Optional]
% 3. Calculate centroid of each spot in GE
% 4. Calculate the pairwise distance between synth and GE spots.

clear; clc
tic
MONOCLINIC_NITI_LATTICE_PARAMS = [2.889 4.12 4.622 90 96.8 90];
CUBIC_NITI_LATTICE_PARAMS      = [3.015 3.015 3.015 90 90 90];
TETRAGONAL_NITI_LATTICE_PARAMS = [3.015 3.015*sqrt(2) 3.015*sqrt(2) 90 90 90];

%% Inputs (distances in mm, angles in degrees)
PHASE                        = 'monoclinic';                               % 'cubic', 'monoclinic' etc.
OPTIONS                      = [0 0 1 0 0];                                % [Lorentz  Rings  Moasic  OrientationDistribution  Frames]
SUPERIMPOSE_MONO_GE          = true;                                       % Superimpose 1st ring from mono GE2 files
SUPERIMPOSE_CUBIC_GE         = false;                                      % Superimpose 1st ring from cubic GE2 files
RING_NUMBER_CUBIC_GE         = 1;                                          % Which cubic ring to superimpose
RING_NUMBER_MONO_GE          = 1;                                          % Which cubic ring to superimpose
IMAGESC_TYPE                 = 'raw';                                      % 'raw': only GE2 data, 'superimposed': GE2 + synth data
MAKE_SYNTH_SCATTER_PLOT      = true;                                       % Plot synth eta, omega pais as scatter plot
MAX_GE_INTENSITY             = 20;                                         % Max intensity of synth spots
PRINT_IMAGES                 = false;                                      % Print GE + synth data to PNGs
CALCULATE_SPOT_MATCH         = false;                                      % Calculate pairwise distance between GE/synth spots
%
MATERIAL_PROPS               = struct();                                   % Material properties (define below)
MATERIAL_PROPS.latticeParams = MONOCLINIC_NITI_LATTICE_PARAMS;                  % [a b c alpha beta gamma] (ang in degrees)
MATERIAL_PROPS.spaceGroup    = 11;
MATERIAL_PROPS.beamEnergy    = 55.318;                                     % In keV
%MATERIAL_PROPS.orientation   = [0.999 1E-5 1E-5 1E-5];
MATERIAL_PROPS.orientation   = [9.925539406453766089e-01        -1.071975186477794745e-01       4.759047213948603382e-03        5.764302537782723529e-02;
                                9.919089671787508777e-01        -1.015517987195288102e-01       1.231356456392354538e-03        7.617293987033339764e-02]; % (June 0 deg) Cubic orientation of the grain (quaternion)
%MATERIAL_PROPS.orientation   = [9.031931474684801175e-01        2.398647726362172772e-01        -3.127018349637547545e-01       -1.700723129273450929e-01]; % (Dec 20 deg) Cubic orientation of the grain (quaternion)
%MATERIAL_PROPS.orientation   = [9.329181291088404215e-01        1.131258077512668075e-01        -3.323692948023018179e-01       7.998104696304578209e-02]; % (Dec 0 deg) Cubic orientation of the grain (quaternion)
MATERIAL_PROPS.variantsLoop  = [1:24];                                       % Variants in the cubic point group to loop over. e.g. [1:24], [10 14]
%
DETECTOR_PROPS               = struct();                                   % Detector params (define below)
DETECTOR_PROPS.distance      = 1110.0;                                 % Detector to sample distance (mm)
%
GE2_PROPS                    = struct();                                   % Properties for read/writing GE2 files
GE2_PROPS.deltaOmega         = 0.1;                                        % Omega step for GE2 files
GE2_PROPS.outPath            = '/home/software-dev/niti_spot_analysis/RCV12_synth';
GE2_PROPS.outName            = 'NiTi_B19_RsynthAB_RCV12_01.ge2';
%%
% Distance from sample to detector (m)
distance = DETECTOR_PROPS.distance * 1e-3;

% delta omega
deltaOmega = GE2_PROPS.deltaOmega;  %(deg)
% deltaOmega = 360;

% Incoming beam wavelength (A)
beamEnergy = MATERIAL_PROPS.beamEnergy;        % (keV)
h = 6.626E-34;
c = 3.000E+08;
e = 1.602E-19;
lambda = h * c / ( 1000 * beamEnergy * e ) * 1e10; % (A)

% .ge2 destination
outpath = GE2_PROPS.outPath;
outname = GE2_PROPS.outName;
l1=[1;0;0]; l2=[0;1;0]; l3=[0;0;1];

a0 = 3.015;
a = MATERIAL_PROPS.latticeParams(1);
b = MATERIAL_PROPS.latticeParams(2);
c = MATERIAL_PROPS.latticeParams(3);
alpha = MATERIAL_PROPS.latticeParams(4) * pi/180.0;
beta = MATERIAL_PROPS.latticeParams(5) * pi/180.0;
gamma = MATERIAL_PROPS.latticeParams(6) * pi/180.0;

% Get matrix representation of 24 cubic symmetry ops
CUBICROTS = getSymmetryOp();

if(SUPERIMPOSE_CUBIC_GE)
    load(['ge_cubic_eta_omega_ring_' num2str(RING_NUMBER_CUBIC_GE)]);
    MAX_GE_INTENSITY = max(ge_cubic_eta_ome_ring_small(:));
    ge_data = ge_cubic_eta_ome_ring_small;
end

if(SUPERIMPOSE_MONO_GE)
    load(['ge_mono_eta_omega_ring_' num2str(RING_NUMBER_MONO_GE)]);
    eval(['ge_mono_eta_ome_ring_small = ge_mono_eta_omega_ring_' num2str(RING_NUMBER_MONO_GE) ';']);
    ge_mono_eta_ome_ring_small = [ge_mono_eta_ome_ring_small; zeros(15, 3600)];
    MAX_GE_INTENSITY = max(ge_mono_eta_ome_ring_small(:));
    ge_data = ge_mono_eta_ome_ring_small;
end

% Calculate centroid of spots in the GE data
if(CALCULATE_SPOT_MATCH)
    ge_data_filtered = imfilter(100*ge_data, fspecial('gaussian', 12, 6));
    ge_data_peaks = FastPeakFind(ge_data_filtered, 3);
    ge_data_peaks_y = ge_data_peaks(1:2:end);
    ge_data_peaks_x = ge_data_peaks(2:2:end);
    spot_distances = [];
end

for jj = MATERIAL_PROPS.variantsLoop

    disp(['Processing variant ' num2str(jj)])
    
    %% List of reflections
    if(SUPERIMPOSE_MONO_GE)
    initial_hkl_list = [0 0 1;
        0 1 1;
        1 0 0;
        -1 0 1;
        1 1 0;
        1 0 1;
        0 0 2;
        -1 1 1;
        0 2 0;
        1 1 1];
    
    initial_hkl_list = initial_hkl_list(RING_NUMBER_MONO_GE, :);
    elseif(SUPERIMPOSE_CUBIC_GE)
        initial_hkl_list = [
        0 1 0;
        1 1 0;
        1 1 1;
        0 2 0;
        1 2 0];
    
    initial_hkl_list = initial_hkl_list(RING_NUMBER_CUBIC_GE, :);
    end
    
    % Get all symmetry ralated Miller indices for a plane
    if(SUPERIMPOSE_MONO_GE)
        temp = [];
        for j = 1:size(initial_hkl_list, 1)
            %vecin = cubic_symmetries( transpose( initial_hkl_list(j,:) ) );
            %vecin = [0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
            vecin = [0 0 1; 0 0 -1];
            temp = vertcat(temp, vecin);
        end
        hkl_list=unique(temp,'rows'); % Mono symmetry related planes
    elseif(SUPERIMPOSE_CUBIC_GE)
        temp = [];
        for j = 1:size(initial_hkl_list, 1)
            % Need to figure out what's wrong with mono_symmetries function
            vecin = cubic_symmetries( transpose( initial_hkl_list(j,:) ) );
            temp = vertcat(temp, vecin);
        end
        hkl_list=unique(temp,'rows'); % Cubic symmetry related planes
    end

    %% Create mosaic geometry - Parallelopiped
    if OPTIONS(3) == 1
        SpecimenSize.x = 1e-3;  % (m)
        SpecimenSize.y = 1e-3; % vertical, beam height???
        SpecimenSize.z = 1e-3; % beam direction
        
        NumBlocks.x = 5;  Blocksize.x = SpecimenSize.x / NumBlocks.x;
        NumBlocks.y = NumBlocks.x;  Blocksize.y = SpecimenSize.y / NumBlocks.y;
        NumBlocks.z = NumBlocks.x;  Blocksize.z = SpecimenSize.z / NumBlocks.z;
        
        Block.x = linspace(-SpecimenSize.x/2+Blocksize.x/2,  SpecimenSize.x/2-Blocksize.x/2,  NumBlocks.x)';
        Block.y = linspace(-SpecimenSize.y/2+Blocksize.y/2,  SpecimenSize.y/2-Blocksize.y/2,  NumBlocks.y)';
        Block.z = linspace(-SpecimenSize.z/2+Blocksize.z/2,  SpecimenSize.z/2-Blocksize.z/2,  NumBlocks.z)';
        
        Volume = (SpecimenSize.x / NumBlocks.x) ^ 3; 
        
        Geometry = zeros(NumBlocks.x*NumBlocks.y*NumBlocks.z, 4);  gg = 1;
        for ix = 1 : NumBlocks.x
            for iy = 1 : NumBlocks.y
                for iz = 1 : NumBlocks.z
                    Geometry(gg,:) = [Block.x(ix) Block.y(iy) Block.z(iz) Volume];
                    gg = gg + 1;
                end
            end
        end
    end
    
    % Geometry = [0 0 0 1];
    %% Orientation
        % RFix is an operation belonging to the cubic point group (24 total)
    RFix = squeeze(CUBICROTS(:, :, jj));
    R45 = [1.0000   -0.0000    0.0000; ...
           0.0000    0.7071    0.7071; ...
           0.0000   -0.7071    0.7071];
    % Loop over all cubic grain orientations
    Lights = []; % Synth data s stored in this matrix
    for rr = 1:size(MATERIAL_PROPS.orientation, 1)
        % Cubic grain orientation
        quat = MATERIAL_PROPS.orientation(rr, :);
        if(strcmpi(PHASE, 'cubic'))
            Orientation = ((quat2rot(quat)) * RFix)';
            disp(Orientation)
        elseif(strcmpi(PHASE, 'monoclinic'))
            Orientation = (( R45 * quat2g(quat) ) * RFix)';
        end

        %% Calculate virtual detector data
        fprintf('Processing ');
        for j = 1:size(hkl_list,1)

            fprintf(['(' num2str(hkl_list(j,:)) '), '])

            Lights_temp = MosaicLightUp_B19_V12(hkl_list(j,:), lambda, distance, Orientation, Geometry, MATERIAL_PROPS);

            if OPTIONS(3) == 0 && OPTIONS(4) == 1
                Lights = vertcat(Lights, Lights_temp1, Lights_temp2);
            else
                Lights = vertcat(Lights, Lights_temp);
            end
        end

        fprintf('\n')

        fprintf('Light Up Finished.\n\n')
        % display(Lights)
    end
    
    if(size(Lights, 1) == 0 || size(Lights, 2) == 0)
        continue;
    end
    
    LightsOriginal = Lights;
    R = distance * tan(Lights(:,2)); IntTemp = Lights(:,4);
    
    % Remove anything with zeroish intensity or that hits outside the detector
    AcceptedIndexes1 = find(abs(R) <= 0.2048);
    Lights = Lights(AcceptedIndexes1, :);
    AcceptedIndexes2 = find(abs(Lights(:,4)) >= 10);
    Lights = Lights(AcceptedIndexes2, :);
    Lights(:,4) = 1*MAX_GE_INTENSITY;  % Set the same intensity to all spots
    
    R = distance * tan(Lights(:,2));
    zeta_x = R .* cos(Lights(:,3));
    zeta_y = R .* sin(Lights(:,3));
    
    
    %% Assign Lorentzian Distribution
    radius = 1e-3;
    if OPTIONS(1) == 1
        LorentzNum = 7;
        LorentzRadius = logspace(log10(radius)-2, log10(radius), LorentzNum)';
        LorentzTheta = linspace(0, 2*pi, LorentzNum*2)';
        
        k = 1;
        for im = 1 : length(LorentzRadius)
            for jm = 1 : length(LorentzTheta)
                
                XM(k,1) = LorentzRadius(im) * cos(LorentzTheta(jm));
                YM(k,1) = LorentzRadius(im) * sin(LorentzTheta(jm));
                
                Rcenter = 0;
                Rwidth = radius * 0.75;
                
                Z(k,1) = (1/pi * 1/2*Rwidth / ( (LorentzRadius(im)-Rcenter)^2 + (1/2*Rwidth)^2 )) * 1e-5;
                k = k + 1;
            end
        end
        
        zeta_xNew = zeta_x;
        zeta_yNew = zeta_y;
        clear zeta_x zeta_y
        zeta_x = [];  zeta_y = [];  Int = [];
        for ii = 1 : length(zeta_xNew)
            zeta_xTemp = XM + zeta_xNew(ii);
            zeta_yTemp = YM + zeta_yNew(ii);
            IntTemp = Lights(ii,4) * Z;
            
            zeta_x = vertcat(zeta_x, zeta_xTemp);
            zeta_y = vertcat(zeta_y, zeta_yTemp);
            Int = vertcat(Int, IntTemp);
        end
        clear LorentzRadius LorentzTheta XM YM Z zeta_xNew zeta_yNew
        
    else
        Int = Lights(:,4);
    end
    
    
    %%
    
    %% Rings
    if OPTIONS(2) == 1
        for i = 1 : size(hkl_list, 1)
            hkl = hkl_list(i,:);  h = hkl(1);  k = hkl(2);  l = hkl(3);
            
            dhkl = 1 / sqrt( h^2/(a^2 * sin(beta)^2) + k^2/b^2 + l^2/(c^2*sin(beta)^2) - 2 * h * l * cos(beta)/(a * c * sin(beta)^2) );
            RingTheta = asin(lambda / (2 * dhkl));  RingTTheta = 2 * RingTheta;
            
            RingRadius = distance * tan(RingTTheta);
            RingAngles = 0 : 0.1 : 2*pi;
            RingX = RingRadius * cos(RingAngles);
            RingY = RingRadius * sin(RingAngles);
            figure(1)
            plot(RingX, RingY, 'r'); hold on
            text(RingRadius+0.00, 0.0, num2str(hkl), 'Rotation', 270, 'Color', 'r'); hold on
            clear RingAngles RingX RingY
        end
    end
    
    
    %% Sort into frames and export tiffs
    if OPTIONS(5) == 1
        h = tic;
        Lights(:,1) = Lights(:,1) * 180/pi + 180;
        
        Omega = 0 : deltaOmega : 360;
        for i = 1 : length(Omega)-1
            DeltaOmegaIndexes = find(Lights(:,1) >= Omega(i)  &  Lights(:,1) < Omega(i+1));
            FrameLights{i,1} = Lights(DeltaOmegaIndexes, :);
            FLights = FrameLights{i,1};
            
            R = distance * tan(FLights(:,2));
            zeta_x = R .* cos(FLights(:,3));
            zeta_y = R .* sin(FLights(:,3));
            Int = FLights(:,4);
            binsize = 200e-6;       % 200 um
            xbins = -0.2048+binsize/2 : binsize : 0.2048-binsize/2;
            ybins = xbins;
            [nx, idx] = histc(zeta_x, xbins);
            [ny, idy] = histc(zeta_y, ybins);
            out = zeros(length(xbins));
            for ii = 1:length(idx)
                %         out(idx(ii),idy(ii)) = out(idx(ii),idy(ii)) + Int(ii); % sum %%%%%%%%%%%%%%%
                out(idx(ii),idy(ii)) = max(out(idx(ii),idy(ii)),Int(ii));
            end
            
            if i == 1
                %         fo = ExportFrame(FrameLights{i,1}, distance, i, 1);
                fo = fopen(fullfile(outpath,outname),'w');
                header = zeros(1,8192);
                fwrite(fo,header,'uint8');
                fseek(fo,size(header,2),'bof');
                fclose(fo);
            end
            %         fo = ExportFrame(FrameLights{i,1}, distance, i, 0);
            fo = fopen(fullfile(outpath,outname),'a+');
            fwrite(fo, out, 'uint16');
            fclose(fo);
            
        end
        disp('GE2 CREATION COMPLETE');
        toc(h)
        
        load handel
        sound(y,Fs)
    end
    
    % Get unique eta, omega, theta solutions (with 0.001 fuzz)
    Lights = uniquetol(Lights, 0.001, 'ByRows', true);
    % MosaicLightUp has angles in [-180 180]. So add 180.
    omega_synth = Lights(:, 1)*180.0/pi + 180.0;
    twotheta_synth = Lights(:, 2);
    eta_synth = Lights(:, 3)*180.0/pi + 180.0;
    
    % Generate an eta-omega map from the (eta, omega, two-theta) list
    etaOmegaMap = zeros(round(360*1.0/GE2_PROPS.deltaOmega));
    % Note the order here - omega: rows, eta: columns
    for kk = 1:numel(omega_synth)
        etaOmegaMap(max(1, round(1.0/GE2_PROPS.deltaOmega*omega_synth(kk))), max(1, round(1.0/GE2_PROPS.deltaOmega*eta_synth(kk)))) = MAX_GE_INTENSITY;
    end
    %etaOmegaMap = etaOmegaMap';
    % Superimpose GE2 data for 1st cubic ring
    if(SUPERIMPOSE_CUBIC_GE)
        etaOmegaMap = etaOmegaMap + ge_cubic_eta_ome_ring_small;
    end
    % Superimpose GE2 data for 1st monoclinic ring
    if(SUPERIMPOSE_MONO_GE)
        etaOmegaMap = etaOmegaMap + ge_mono_eta_ome_ring_small;
    end
    %%
    % Plotting
    h(jj) = figure;
    
    if(strcmpi(IMAGESC_TYPE, 'raw'))
        % Plot only the raw GE2 data and not the synth data in imagesc plot
        if(SUPERIMPOSE_CUBIC_GE)
            imagesc((ge_cubic_eta_ome_ring_small)); axis equal;
        elseif(SUPERIMPOSE_MONO_GE)
            imagesc((ge_mono_eta_ome_ring_small)); axis equal;
        else
        end
           
    else
        % Plot raw data + synth data (superimposed)
        imagesc(etaOmegaMap); axis equal;
    end
    
    hold on;
    % Plot synth eta, omega as a scatter plot
    if(MAKE_SYNTH_SCATTER_PLOT)
        scatter(1.0/GE2_PROPS.deltaOmega*eta_synth, 1.0/GE2_PROPS.deltaOmega*omega_synth, 16, 'r');
    end
    % Print variant number
    text(10, 10, num2str(jj), 'Color', 'w', 'FontSize', 20);
    hold off;
    drawnow;
    if(PRINT_IMAGES)
        print(h(jj), ['plots/synth_verification_mono_CV' num2str(jj)], '-dpng');
    end
    
    if(CALCULATE_SPOT_MATCH)
        debug_this_block = false;
        if(debug_this_block)
            figure;
            hold on;
            scatter(ge_data_peaks_y, 360 - ge_data_peaks_x, 36, 'b', 'filled');
            scatter(eta_synth, 360 - omega_synth, 36, 'r', 'filled');
            hold off;
            xlim([0 360]); ylim([0 360]);
            axis square;
            drawnow;
        end
        %
        eta_ome_synth_data = [eta_synth (360 - omega_synth)];
        eta_ome_ge_data = [ge_data_peaks_y (360 - ge_data_peaks_x)];
        eta_ome_spot_distances = pdist2(eta_ome_ge_data, eta_ome_synth_data);
        eta_ome_spot_distances_min = min(eta_ome_spot_distances, [], 2);
        spot_distances = [spot_distances; eta_ome_spot_distances_min'];
    end
end

clearvars -except a b c beta a0 RFix Orientation omega_synth twotheta_synth ...
                  eta_synth jj CUBICROTS ge_cubic_eta_ome_ring_small ...
                  spot_distances Geometry

toc