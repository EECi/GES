function [Nu, varargout] = lamorturb(Gr, Re, varargin)

% This function calculates the Nusselt number due to free and forced
% convection, and then chooses the larger of the two. In each case the
% convection may be laminar or turbulent, and this depends on the magnitude
% of the Grashof and Reynolds number respectively.
%
% The dimensionless number Sh, used in computing the heat exchange due to
% evaporation and condensastion of water, depends on the Lewis number of
% air in a manner that depends on whether convection is free or forced. It
% is treated as an optional argument.

% Free convection, laminar or turbulent
free = Gr < 1e5;
Nu_G = 0.5*free.*Gr.^0.25 + 0.13*(1 - free).*Gr.^0.33;

% Forced convection, laminar or turbulent
forced = Re < 2e4;
Nu_R = 0.6*forced.*Re.^0.5 + 0.032*(1 - forced).*Re.^0.8;

% Calculating the larger of the two, and Sh
x = Nu_G > Nu_R;
Nu = x.*Nu_G + (1 - x).*Nu_R;

if nargin == 3
    Le = varargin{1};
    Sh = x.*Nu*Le^0.25 + (1 - x).*Nu*Le^0.33;
    varargout{1} = Sh;
end
    
end