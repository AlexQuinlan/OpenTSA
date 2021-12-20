function dToT_mag = Eval_dToTZS(T0, Tcomp, Ex, Ex_0, Zstar, Kx, delS_LC)
%% Zhang and Sandor model
    %Input arguments
    % T0 - mean temperature of sample
    % Tcomp - temperature compensation toggle
    % Ex - Current Laminate longitudinal stiffness
    % Ex_0 - Undamaged laminate longitudinal stiffness
    % Zstar - adjustment factor
    % Kx - thermoelastic constant for the laminate
    % delS_LC - Stress amplitude (in MPa)
    
    
    if ~Tcomp  % if no temperature compensation, set terms to 1
        Cstar = 1;
    else
        Cstar = CpT0(300) / CpT0 (T0);  
    end
    
    
    D = (Ex_0 - Ex) / Ex_0;
    dToT_mag = Zstar * Cstar * Kx * delS_LC / (1-D)^2;
    
end