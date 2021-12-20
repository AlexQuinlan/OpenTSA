function [tauPly_cmplx,tauSurf_cmplx] = Burgers_1D_engine(f, gamma, surfthick, ...
    ply_thickness, hp, xdisc, t_disc, cycles)
%%%%%%%%% Input Arguments %%%%%%%%%%%%%
%%%%% Required
% f               - frequency in Hz
% gamma           - thermal diffusivity
% surfthick       - Thickness of the SRRL
% ply_thickness   - Thickness of each ply
% hp              - 1/2 of the ply count (assuming symmetry)
%%%%% Optional
% xdisc           - Discretization of thickness from surface to mid-plane
% t_disc          - Discretization of time in each per period  
% cycles          - Number of cycles until steady state achieved
%%%%% Default Settings
if ~exist('xdisc','var')
  xdisc = 481;
end
if ~exist('t_disc','var')
  t_disc = 100;
end
if ~exist('cycles','var')
  cycles = 10;
end

%%%%%%%%%% Initialization %%%%%%%%%%%%%%
pf = @(y) y;  ff = @(y) y;  % For sine fitting function
m=0;  % For solver (symmetry parameter?)
x = linspace(0,surfthick+ply_thickness*hp,xdisc); 
w=f*2*pi;  % radians
t = linspace(0,cycles/(f),t_disc*cycles);
refsig = cos(w*t);
tausMag_fp=zeros(hp,1);
tausPh_fp=zeros(hp,1);
%%%%%%%%%% Start Solver %%%%%%%%%%
%tic   % Start time for optimization
for i_ply = 0:hp  % 
    fnh = @(x,t,u,d) heatcyl(x,t,u,d,w,gamma,i_ply,ply_thickness,surfthick);  
    sol = pdepe(m,fnh,@heatic,@heatbc,x,t);
    tauSurf = sol(:,1,1);  % tau at the surface
    dtau = range(tauSurf)/2;
    [~,ph,~,~]=SineFittingLeastSquare(f,t, tauSurf ,pf,ff);
    if i_ply == 0
        tausMag_surf = dtau;
        tausPh_surf  = ph;
    else
        tausMag_fp(i_ply) = dtau;
        tausPh_fp(i_ply) = ph;
    end
end  % End ply loop

tauPly_cmplx = tausMag_fp .* exp(1i * tausPh_fp);
tauSurf_cmplx = tausMag_surf .* exp(1i * tausPh_surf);

end  %  End function
%%%%%%%%%  Functions %%%%%%%%%%%%
% Functions are called from the pdepe solver, and cannot recieve any
% additional arguments.  Thus 'evalin' is used to obtain variables
% from the matlab workspace
function [c,f,s] = heatcyl(x,t,u,dudx,w,gamma,i_ply,ply_thickness,surfthick)
%      w=evalin('base','w');   
%      gamma=evalin('base','gamma');
%     load('temp_func_var.mat')
    c = 1;
    f = dudx*gamma;
    Hj = (1+floor((x-surfthick)/ply_thickness)) == i_ply;
    s=Hj*cos(w*t)*w;
end

function u0 = heatic(x)
    % Set Initial Conditions
    u0=0;
end

function [pL,qL,pR,qR] = heatbc(xl,ul,xr,ur,t)
    % Set Boundary Conditions
    pL = 0; % value of dtau/dz at left
    qL = 1; % multiply by dtau/dz
    pR = 0; % value of dtau/dz at right
    qR = 1; 
end

function Hj = funcHj(x)
%     load('temp_func_var.mat')
     curPly=evalin('base','i_ply');
     plythick=evalin('base','ply_thickness');
     surfthick=evalin('base','surfthick');
%     curPly = i_ply;
%     plythick = ply_thickness;
    
    if (1+floor((x-surfthick)/plythick)) == curPly
        Hj = 1;
    else
        Hj = 0;
    end
end
