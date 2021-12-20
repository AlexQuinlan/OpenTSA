classdef TSAsample < handle
   properties
       % Experimental
       T0
       delepsx
       delepsy
       delSigx
       f
       exp_dT
       exp_dToT
       cyc
       Ex_exp
       % Laminate info
       ply_thickness
       laminfo
       mpf
       layup
       thickness
       hp
       surfthick
       plies % custom class
       
       % Material Property Groups
       pp  % ply properties
       rp  % resin properties
       ABD % Calculated laminate matrix
      
       % Conduction Contribution Factors
       tauSRRL    
       tauPlies
       gamma
       
       % Other
       pde_params
       
       % Theoretical Predictions
       dToT_eq20  % result from current paper equation 20
       dToT_eq20cmplx % Complex result from eq20
       dToT_ZS    % result from Zhang & Sandor equation
       dToT_ZSmod % result from modified Zhang & Sandor
       dToT_SRRL  %
       dToT_GLOB
       dToT_SURFP
       
       dmg1 % ply-by-ply vector of damage in direction 1
       dmg2 % ply-by-ply vector of damage in direction 2
       DMG  % laminate scalar determined by stiffness loss
   end
   
   methods
       
       %initialize
       function obj=TSAsample(laminfo, T0, dex, dey, dSx, f, surfthick )
           
           %laminfo - laminate information file
           obj.laminfo=laminfo;
           %T0 - experimental mean temperature
           obj.T0 = T0;
           %delepsx - global x strain range
           obj.delepsx = dex;
           %delepsy - global y strain range
           obj.delepsy = dey;
           %delSigx - global x stress range
           obj.delSigx = dSx;
           % f - loading frequency
           obj.f = f;
           % surfthick - surface thickness
           obj.surfthick  = surfthick;
           
           % default parameters - change manually if desired
           obj.pde_params = [481, 100, 15]; % [xdisc, t_disc, cycles]
           
           
       end
       
       function build_plies(obj, dmg1, dmg2)  
           
          if exist('dmg1','var') 
              obj.apply_damage(dmg1,dmg2);
          end
          zbot = -1/2*sum(obj.thickness);
          for i_ply = 1:length(obj.layup)
             pplies{i_ply} = AQply(obj.pp.E_1, obj.pp.E_2, obj.pp.G_12, ...
                 obj.pp.nu_12,obj.pp.nu_21);
           %%  plies{i_ply}.calcDmgProps(fiberDmg(i_ply), matrixDmg(i_ply));
             % Need to add in damage later
             pplies{i_ply}.calcDmgProps(obj.dmg1(i_ply), obj.dmg2(i_ply)); % No damage applied
             pplies{i_ply}.makePlyQ();
             pplies{i_ply}.placePly(obj.layup(i_ply), obj.thickness(i_ply), zbot)
             pplies{i_ply}.makeTmatrix();
             pplies{i_ply}.calcQxys();
             zbot = zbot + obj.thickness(i_ply);
          end
          
          obj.plies = pplies;
       end
       
       
       function read_laminfo(obj)
           run([obj.laminfo])  % run code program
           obj.layup = layup;
           obj.thickness = thickness;
           obj.hp = ceil(length(layup)/2);
           obj.ply_thickness = mean(thickness);
           %mpf - material property file
           obj.mpf=material_File;
           obj.dmg1 = zeros(1,length(layup)); % Default to undamaged
           obj.dmg2 = zeros(1,length(layup)); % Default to undamaged
       end
       
       function read_mpf(obj)
           run([obj.mpf])  % run code program
           obj.pp = pp;
           obj.rp = rp;
           obj.gamma = (obj.pp.k) / (obj.pp.rho * obj.pp.Cp);
       end
       
       function calc_tau(obj)
           [tauPly, tauSurf] = Burgers_1D_engine(obj.f, obj.gamma, ...
               obj.surfthick, obj.ply_thickness, obj.hp, ...
               obj.pde_params(1), obj.pde_params(2), obj.pde_params(3) );
           obj.tauPlies = tauPly;
           obj.tauSRRL = tauSurf;
           
       end
       
       function set_experimental(obj, dT)
          obj.exp_dT = dT;
          obj.exp_dToT = dT ./ obj.T0;
       end
       
       function evaluate_dToTAQ(obj)
           [dToT_mag, dToT_cplx] = Eval_dToTAQ(obj.T0, obj.plies, ...
               obj.delepsx, obj.delepsy, obj.pp, obj.rp, obj.tauPlies, ...
               obj.tauSRRL, 1);
            obj.dToT_eq20 = dToT_mag;
            obj.dToT_eq20cmplx = dToT_cplx;
       end
       
       function evaluate_SRRL(obj)
           obj.dToT_SRRL = obj.rp.Kr * (obj.rp.Er/(1-obj.rp.nur) * ...
               (obj.delepsx + obj.delepsy));   
       end
       
       function evaluate_surfply(obj)                    
            p1=obj.plies{1};
            obj.dToT_SURFP = abs( ...
                (p1.Q126(1,1)*obj.pp.K1 + ... 
                p1.Q126(1,2)*obj.pp.K2)*p1.eps126(1) + ...
                (p1.Q126(1,2)*obj.pp.K1 + ...
                p1.Q126(2,2)*obj.pp.K2)*p1.eps126(2)  );
       end
       
       function evaluate_global(obj)
           
           obj.dToT_GLOB = abs(  ...
               [obj.pp.Kx, obj.pp.Ky, 0]*obj.ABD(1:3,1:3)/sum(obj.thickness) ...
                * [obj.delepsx ; obj.delepsy; 0]  ); 
            
       end
       
       function calc_ABD(obj)
           Ap = zeros(3,3,length(obj.plies));
           Bp = zeros(3,3,length(obj.plies));
           Dp = zeros(3,3,length(obj.plies));    
           for i_ply = 1:length(obj.plies)
               p = obj.plies{i_ply};  %shorthand for space
               Ap(:,:,i_ply) = p.Qxys * (p.z2 - p.z1);
               Bp(:,:,i_ply) = 1/2 * p.Qxys * (p.z2^2 - p.z1^2);
               Dp(:,:,i_ply) = 1/3 * p.Qxys * (p.z2^3 - p.z1^3);
           end    
           A = sum(Ap,3); B=sum(Bp,3); D=sum(Dp,3);
           obj.ABD = [A, B; B,D];
       end
       
       function apply_damage(obj, df, dm)
           %% df - vector of fiber damage in each ply
           %% dm - vector of matrix damage in each ply
           
           if length(df) < length(obj.layup)
              obj.dmg1=[df flip(df)];
              obj.dmg2=[dm flip(dm)];
           else
              obj.dmg1=[df ];
              obj.dmg2=[dm ];
           end
%            for i_p = 1:length(obj.plies)
%                 obj.plies{i_p}.df = D1(i_p);
%                 obj.plies{i_p}.dm = D2(i_p);
%            end
          
           
       end
       
       function calc_globalCTE(obj)         
           Ap = zeros(3,3,length(obj.plies));
           Bp = zeros(3,3,length(obj.plies));
           Dp = zeros(3,3,length(obj.plies));
           pabar = zeros(3, length(obj.plies));         
           for i_ply = 1:length(obj.plies)
               pabar(:,i_ply) = obj.plies{i_ply}.t * ( inv(obj.plies{i_ply}.T)...
                   * obj.plies{i_ply}.Q126 * [obj.pp.alpha1 ; obj.pp.alpha2 ; 0]);
               p = obj.plies{i_ply};  %shorthand for space
               Ap(:,:,i_ply) = p.Qxys * (p.z2 - p.z1);
               Bp(:,:,i_ply) = 1/2 * p.Qxys * (p.z2^2 - p.z1^2);
               Dp(:,:,i_ply) = 1/3 * p.Qxys * (p.z2^3 - p.z1^3);
           end
           A = sum(Ap,3); B=sum(Bp,3); D=sum(Dp,3);
           obj.ABD = [A, B; B,D];
           abd = inv(obj.ABD);            
           alphaGlob = abd(1:3,1:3)*sum(pabar,2);  
           obj.pp.alpha_x = alphaGlob(1);
           obj.pp.alpha_y = alphaGlob(2);
            
       end
       
       function calc_Kxy(obj)           
           obj.pp.Kx = obj.pp.alpha_x / (obj.pp.Cp * obj.pp.rho);
           obj.pp.Ky = obj.pp.alpha_y / (obj.pp.Cp * obj.pp.rho);
       end
       
       function calc_globprops(obj)
          abd = pinv(obj.ABD);
          obj.pp.nuxy = -abd(2,1)/abd(1,1);
          obj.pp.Ex = 1/(sum(obj.thickness)) *...
              (obj.ABD(1,1)-(obj.ABD(1,2)^2/obj.ABD(2,2)));
          obj.pp.Ey = 1/(sum(obj.thickness)) *...
              (obj.ABD(2,2)-(obj.ABD(1,2)^2/obj.ABD(1,1)));
       end
       
       function evaluate_ZS(obj, Ex_0)
           % Without Temperature Compensation
           obj.dToT_ZS = Eval_dToTZS(obj.T0, 0, obj.Ex_exp, Ex_0,...
               1, obj.pp.Kx, obj.delSigx);
           % With Temperature Compensation
           obj.dToT_ZSmod=Eval_dToTZS(obj.T0, 1, obj.Ex_exp, Ex_0,...
               1, obj.pp.Kx, obj.delSigx);
       end
       
       function plotLam(obj, figno, x, y, exper)
           figure(figno)
           
           %% x-axis
           switch x
               case 'freq'
                    xval = obj.f;
                    xlabel('Frequency (Hz)')
               case 'cyc'
                   xval = obj.cyc;
                   xlabel('Cycles')
               otherwise
                   disp("Invalid X case")
           end
           
           switch y
               case 'dT'
                    expval = obj.exp_dT;
                    ylabel('$\Delta T (mK)$')
               case 'dToT'
                    expval = obj.exp_dToT;
                    e20val = obj.dToT_eq20;
                    surfval = obj.dToT_SURFP;
                    srrlval = obj.dToT_SRRL;
                    globval = obj.dToT_GLOB;
                    ZSval = obj.dToT_ZS;
                    ZSmodval = obj.dToT_ZSmod;
                    ylabel('$\Delta T/T (K/K)$')
               otherwise
                    disp("Invalid X case")
           end
           
           plot(xval , expval, 'r*')
           hold on
           plot(xval , e20val, 'ok')
           plot(xval , obj.dToT_SURFP, 'sb')
           plot(xval , obj.dToT_SRRL, '^b')
           plot(xval , obj.dToT_GLOB, 'vb')
           plot(xval , obj.dToT_ZS, '^g')
           plot(xval , obj.dToT_ZSmod, 'vg')
           
           legend('Experimental','eqn. (20)','Surface Ply',...
               'SRRL', 'Global','ZS','ZS mod')
       end
       
       
   end
end
       
       