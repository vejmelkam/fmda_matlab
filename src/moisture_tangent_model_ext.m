

%
% Tangent linear model to the moisture_model_ext.m extended nonlinear 
% moisture model.
%
%  This function will run the tangent linear moisture model for one time step.
%
%  Synopsis: 
%            m_ext = moisture_model_ext_tangent(T, Tk, Q, P, m_ext, f_info, w, r, dt)
%
%  Arguments:
%
%            T - the temperature in Kelvin
%            Tk - the nominal time lags for each fuel class (seconds)
%            Q - the water vapor content (fraction, dimensionless)
%            P - the current surface atmospheric pressure (Pascals)
%            m - the extended state (see moisture_model_ext.m) including
%                current moisture content in the fuel (dimless. fraction)
%            f_info - the fuel type and location index for each fuel
%                   modeled in m_ext (1st col = f_type, 2nd col = f_loc)
%            r - the current rainfall intensity (mm/h, valid for current
%                          time step)
%            dt - the integration step [s]
%
%
%  Returns:
%
%            Jm_ext - the Jacobian of the nonlinear model at previous (m_ext)
%

function Jm_ext = moisture_tangent_model_ext(Tk, Ew, Ed, m_ext, r, dt)

    k = length(Tk);                 % number of fuel classes
    r0 = 0.05;                      % threshold rainfall [mm/h]
    rk = 8;                         % saturation rain intensity [mm/h]
    Trk = 14 * 3600;                % time constant for wetting model [s]
    
    % first, we break the state vector into components
    m = m_ext(1:k);
    dlt_E = m_ext(k+1);
    
    % rescale to fractions
    % modification: extended model equilibria affected by assimilation
    Ed = Ed + dlt_E;
    Ew = Ew + dlt_E;
    
    % saturation model, equi and rlag are specific to fuel type and
    % location
    rlags = zeros(k+2,1);
    model_ids = zeros(k,1);
    
    % if rain model is active
    if(r > r0)
        model_ids(:) = 3;
        rlags(1:k) = 1.0 / Trk .* (1 - exp(-(r - r0) / rk));
    else
        model_ids(:) = 4;
        rlags(1:k) = 1.0 ./ (Tk * 3600);
        
        model_ids(m < Ew) = 2;
        model_ids(m > Ed) = 1;
    end
    
    % select appropriate integration method (small change -> big errors in
    % the standard solution)
    change = dt * rlags;
    Jm_ext = zeros(k+2);
    for i=1:k
        
        if(change(i) >= 0.01)
            
            % partial m_i/partial m_i
            if(model_ids(i) < 4)
                Jm_ext(i,i) = exp(-change(i));
            else
                Jm_ext(i,i) = 1.0;
            end
            
            % precompute partial m_i/partial equi
            dmi_dequi = (1.0 - exp(-change(i)));

        else
            
            % partial dm_i/partial m_i
            if(model_ids(i) < 4)
                Jm_ext(i,i) = 1.0 - change(i) * (1 - 0.5 * change(i));
            else
                Jm_ext(i,i) = 1.0;
            end
            
            % partial m_i/partial equi
            dmi_dequi = change(i) * (1 - 0.5 * change(i));
                        
        end
        
        % branch according to the currently active model
        if(r <= r0)
            
            % drying/wetting model active
            if((m(i) > Ed) || (m(i) < Ew))

                % partial m_i/partial delta_E
                Jm_ext(i,k+1) = dmi_dequi;
            end

        else
            % rain model active

            % partial m_i/partial deltaS
            Jm_ext(i,k+2) = dmi_dequi;

        end
                
    end
    
    % the equilibrium constants 
    Jm_ext(k+1, k+1) = 1.0;
    Jm_ext(k+2, k+2) = 1.0;
