function Lights = MosaicLightUp_B19_V12(hkl, lambda, distance, R_C2S, Geometry, MATERIAL_PROPS)
% Geometry is a 4-column (x,y,z,V) matrix of center of mosaic blocks
I=eye(3); e1=[1;0;0]; e2=[0;1;0]; e3=[0;0;1];


%% Lab frame
Lx = [1; 0; 0];
Ly = [0; 1; 0];
Lz = [0; 0; 1];

ki_L = 2*pi / lambda * (-Lz);                                               % Incoming beam wave vector in the lab frame

%% Crystal frame
% lattice params
a0 = 3.015;
a = MATERIAL_PROPS.latticeParams(1);
b = MATERIAL_PROPS.latticeParams(2);
c = MATERIAL_PROPS.latticeParams(3);
alpha = MATERIAL_PROPS.latticeParams(4) * pi/180.0;
beta = MATERIAL_PROPS.latticeParams(5) * pi/180.0;
gamma = MATERIAL_PROPS.latticeParams(6) * pi/180.0;

% basis vectors (direct lattice frame)
xtal1 = a * [1; 0; 0];
xtal2 = b * [0; 1; 0];
xtal3 = c * [cos(beta); 0; sin(beta)];

% reciprocal lattice basis vectors  (2*PI FACTOR?)  *********
volume = dot(xtal1, cross(xtal2, xtal3));
rec1 = cross(xtal2, xtal3) / volume *2*pi;
rec2 = cross(xtal3, xtal1) / volume *2*pi;
rec3 = cross(xtal1, xtal2) / volume *2*pi;

%% hkl
ghkl_C = hkl(1)*rec1 + hkl(2)*rec2 + hkl(3)*rec3;                           % reciprocal lattice vector in the xtal frame
ghkl_S = R_C2S * ghkl_C;                                                    % reciprocal lattice vector in the sample frame

%% omegas
% The roots (omega) of this equ are used to build a list of diffraction
% events (omega positions) from all reciprocal-lattice vectors.
syms w
gx_S = dot(ghkl_S, Lx);
gz_S = dot(ghkl_S, Lz);
Equ18 = cos(w)* gz_S - sin(w)*gx_S - lambda*norm(ghkl_S)^2 /4/pi;           % Eq. 18 of Pagan/Miller
omega = eval(solve(Equ18, w, 'Real', true));                                % IMPROVE FOR TIME???

Lights = zeros(2, 4);
zeta0 = zeros(2, 2);
zeta = [];
Int_temp = 0;
Int = 0;
OmegaAll = [];
if isempty(omega) == 0
    ko_L = cell(2, 1);
    
    for ii = 1:2
        OmegaTemp = zeros(size(Geometry, 1), 1) + omega(ii);
        OmegaAll = vertcat(OmegaAll, OmegaTemp);
        
        R_S2L = [cos(omega(ii)) 0 sin(omega(ii));
            0 1 0;
            -sin(omega(ii)) 0 cos(omega(ii))];                              % to sample frame to lab frame
        
        ghkl_L = R_S2L * ghkl_S;
        
        ko_L{ii, 1} = ki_L + R_S2L*ghkl_S;                                  % outgoing/scattered beam vector
        
        ko_L_unit = unit(ko_L{ii, 1});                                      % magnitude
        kox_L = dot(ko_L_unit, Lx);
        koy_L = dot(ko_L_unit, Ly);
        koz_L = dot(ko_L_unit, Lz);
        
        %% Intensity
        flux = 5.5e9;                                                       % (200mA) (p/s)
        r0 = 2.81794e-5;                                                    % (A)
        volume_prime = dot(rec1, cross(rec2, rec3));
        
        % ignoring some stuff
        theta2_Bragg = ang(ko_L{ii, 1}, ki_L);
        UCSF_squared = unitcellstructurefactor_mono(hkl, theta2_Bragg(1,1)/2, lambda);
        N = 1;% N = V / det([xtal1 xtal2 xtal3]);   % number of unit cells per mosaic module
        Int_temp = flux * r0^2 * N * volume_prime * UCSF_squared;  %%%%%%%%%%%%%%%
        
        %%
        zeta0(ii,:) = [-(kox_L/koz_L)*distance  -(koy_L/koz_L)*distance]; % Pagan/Miller Eq. 10
        
        for gg = 1 : size(Geometry, 1)
            p = R_S2L * R_C2S * [Geometry(gg,1); Geometry(gg,2); Geometry(gg,3)];
            px = dot(p, Lx);
            py = dot(p, Ly);
            pz = dot(p, Lz);
            zeta_temp = zeta0(ii,:) + [px-(kox_L/koz_L)*pz   py-(koy_L/koz_L)*pz]; 
            
            XM = zeta_temp(1,1);  YM = zeta_temp(1,2);
            zeta = vertcat(zeta, [XM YM]);
            Int = vertcat(Int, Int_temp);
        end
%         zeta(ii,:) = zeta0(ii,:);
    end
    R = (zeta(:,1).^2 + zeta(:,2).^2).^0.5;
    tth = atan2(R, distance);
    eta = atan2(zeta(:,2), zeta(:,1));
    
    Int = zeros(size(zeta,1), 1) + Int_temp;
    Lights = [OmegaAll tth eta Int];
else
    Lights = [];
end