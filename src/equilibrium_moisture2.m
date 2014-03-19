
%
%  Compute the equilibrium moisture field from station data, that is
%  relative humidity and air temperature.
%
%  synopsis: [Ed,Ew] = equilibrium_moisture(Rh, T)
%
%                      Rh - relative humidity [-]
%                      T  - surface temperature [K]
%                      Ed, Ew - drying/wetting equilibrium [dimensionless,
%                                                fraction between 0 and 1]
%
%  Note: shape of Ed, Ew will be the same as the shape of Rh,T.  The
%  shapes of inputs Rh,T must be equal.
%
%

function [Ed,Ew] = equilibrium_moisture2(Rh, T)

    % relative humidity (percent, at each location, size n x 1)
    H = 100 * Rh;
    
    % drying/wetting fuel equilibrium moisture contents (location specific,
    % n x 1)
    Ed = 0.924*H.^0.679 + 0.000499*exp(0.1*H) + 0.18*(21.1 + 273.15 - T).*(1 - exp(-0.115*H));
    Ew = 0.618*H.^0.753 + 0.000454*exp(0.1*H) + 0.18*(21.1 + 273.15 - T).*(1 - exp(-0.115*H));

    % remap values
    Ed = Ed * 0.01;
    Ew = Ew * 0.01;