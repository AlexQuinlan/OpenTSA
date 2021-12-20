function [dToT_mag, dToT_cplx] = Eval_dToTAQ(T0, plies, del_epsx, del_epsy, ...
                    pp, rp, normtau_cmplx, tausrrl, Astar)

        % Input Arguments
        % T0 - mean temperature during tsa reading
        % plies - structure holding stiffness matrix and strains
        % del_epsx - global longitudinal strain Range
        % pp - ply properties
        % rp - resin properties
        % normtau_cmplx - ply complex tau values**
        % tausrrl - srrl complex tau values**
        % Astar - adjustment factor
        
        % ** In the future, read these in from a specified input file
 
        %% load material properties
        %run(matpropFile)
        K1 = pp.K1;   K2 = pp.K2;   K2f = pp.K2;
        disp("No K2f used.  Add, please");
        Kr = rp.Kr;   Er = rp.Er;   nur = rp.nur;
        
        % Evaluate material property compensation terms
        disp("NOTICE: ACp and Ba2 temporarily disabled")
        ACp=1;Ba2=1;
%         ACp =  CpT0(300) / CpT0 (T0);                
%         Ba2 = cte2T0(T0) / cte2T0(300);
%         
        % Ba1 assumed to be 1 for the CASMaT GFRP material
                
        plydToT = {};  % structure for ply response data
        % Expected thermoelastic heating in Ply 1 (0 degrees)
        plydToT{1} = ACp*(...
            (plies{1}.Q126(1,1)*K1 + plies{1}.Q126(1,2)*K2f*Ba2) * del_epsx + ...
            (plies{1}.Q126(1,2)*K1 + plies{1}.Q126(2,2)*K2f*Ba2) * del_epsy );

        % Expected thermoelastic response in Ply 2 (60 degrees) 
        eps_p2 = plies{2}.T' * [del_epsx; del_epsy; 0];   % get the local ply strains
        
        plydToT{2} = ACp*(...
            (plies{2}.Q126(1,1)*K1 + plies{2}.Q126(1,2)*K2*Ba2) * eps_p2(1) + ...
            (plies{2}.Q126(1,2)*K1 + plies{2}.Q126(2,2)*K2*Ba2) * eps_p2(2) );

        % Expected thermoelastic response in the SRRL (isotropic)
        srrldToT = ACp*Ba2*Kr*(...
            (del_epsx + del_epsy) * (Er / (1-nur)) );
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Can we loop through the plies
        hp = ceil( length(plies) /2 );
        plydToT=zeros(1,length(hp)); % Vector for ply response data
        for i_p = 1:hp       
            if i_p == 1
                K2 = pp.K2;  % change to K2f later
            else
                K2 = pp.K2;
            end
            plies{i_p}.eps126 = plies{i_p}.T' * [del_epsx; del_epsy; 0];   % get the local ply strains
            plydToT(i_p) = ACp*(...
                (plies{i_p}.Q126(1,1)*K1 + plies{i_p}.Q126(1,2)*K2*Ba2) * plies{i_p}.eps126(1) + ...
                (plies{i_p}.Q126(1,2)*K1 + plies{i_p}.Q126(2,2)*K2*Ba2) * plies{i_p}.eps126(2) );
        end
        %%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%      
        tau1 = normtau_cmplx(1);
        tau2 = normtau_cmplx(2);
        %tausrrl = tausMag_surf(tauidx)*exp(1i*tausPh_surf(tauidx));
        
%         dToT_cplx = plydToT{1} * tau1 + ...
%                     plydToT{2} * tau2 + ...
%                     srrldToT * tausrrl ;
         
        dToT_cplx = [plydToT, srrldToT]*[normtau_cmplx ; tausrrl];
        
        dToT_mag = Astar*abs(dToT_cplx);
        






end