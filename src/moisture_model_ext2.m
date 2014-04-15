%
% Modified and extended moisture model simulation (based on model in paper
% Kochanski et al., 2012).  Code structure is based on moisture_model.m
%
%
%
%  This function will run the moisture model for one time step.
%  Synopsis: 
%            m_ext = moisture_model_ext2(Tk, Ew, Ed, m_ext, r, dt,S,k,mdE)
%  Arguments:
%
%            Tk - the nominal time lags for each fuel class [hours]
%            Ed - the drying equilibria
%            Ew - the wetting equilibria
%            m_ext - the extended state of the new moisture model, which
%                     includes the fuel moisture (dimensionless fraction),
%                     and equilibria (E, S - all dimensionless fractions).
%                     Total length is k+2, where k is the no. of fuel
%                     classes
%            r - the current rainfall intensity (mm/h, valid for current
%                          time step) [differs between spatial locations,
%                          vector n x 1]
%            dt - the integration step [s]
%
%
%  Returns:
%
%            m_ext - the new moisture values (and state of extended model)
%                    for the time t+dt
%
%            model_id - an integer identifier of the model that is active
%                      (1 - drying, 2 - wetting, 3 - rain, 4 - dead zone)
%

function [m_ext, model_ids] = moisture_model_ext2(Tk, Ed, Ew, m_ext, r, dt,S,rk,r0,Trk,mdE)

    k = length(Tk);                 % number of fuel classes
    Trk = Trk * 3600;               % time constant for wetting model [s]
    
    % first, we break the state vector into components
    m = m_ext(1:k);
    dlt_E = m_ext(k+1);
    dlt_S = m_ext(k+2);
    
    Ed = max(Ed + dlt_E,0);
    Ew = max(Ew + dlt_E,0);
   
    % where rainfall is above threshold (spatially different), apply
    % saturation model, equi and rlag are specific to fuel type and
    % location
    equi = zeros(k+1,1);
    equi(1:k) = m;
    equi(k+1) = mdE;
    rlags = zeros(k+1,1);
    model_ids = zeros(k,1);
    
    % if rain model is active
    if(r > r0)
        model_ids(:) = 3;
        rlags(1:k) = 1.0 / Trk .* (1 - exp(-(r - r0) / rk));
        equi(1:k) = max(S + dlt_S,0);
    else
        model_ids(:) = 4;
        rlags(1:k) = 1.0 ./ (Tk * 3600);
        
        model_ids(m < Ew) = 2;
        equi(m < Ew) = Ew;
        model_ids(m > Ed) = 1;
        equi(m > Ed) = Ed;
    end
    
    % experimental 5hrs time constant for return to mean dE
    rlags(k+1) = 1.0 / (5 * 3600);
    
    % select appropriate integration method (small change -> big errors in
    % the standard solution)
    changes = dt * rlags;
    for i=1:k+1
        if(changes(i) > 0.01)
            m_ext(i) = m_ext(i) + (equi(i) - m_ext(i)) * (1 - exp(-changes(i)));
        else
            m_ext(i) = m_ext(i) + (equi(i) - m_ext(i)) * changes(i) * (1 - 0.5 * changes(i));
        end
    end
    
