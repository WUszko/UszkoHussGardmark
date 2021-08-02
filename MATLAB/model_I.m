%% Community I
% Two unstructured consumer species competing for two resources
% For units and references, see Table S1.2 in Appendix S1
% Created by Wojciech Uszko (2021)

function dYdt = model_I(t, Y)

%% Body masses (ng dry weight)

B_CS = 100;     % small consumer
B_CL = 1000;    % large consumer

%% Temperature- or body mass-independent parameters

deltaRS = 0.1;  % small resource supply rate
deltaRL = 0.1;  % large resource supply rate

q       = 0;    % functional response (Hill) exponent; if =0 then type II

p       = 0.85; % diet preference
pCSRS   = p;
pCSRL   = 1-pCSRS;
pCLRS   = 1-pCSRS;
pCLRL   = pCSRS;

betaCS  = 0.6;  % small consumer conversion efficiency
betaCL  = 0.6;  % large consumer conversion efficiency

HCSRS   = 0.2;  % half-saturation constant
HCSRL   = 0.2;  % half-saturation constant
HCLRS   = 0.2;  % half-saturation constant
HCLRL   = 0.2;  % half-saturation constant

muCS    = 0.01; % small consumer background mortality rate
muCL    = 0.01; % large consumer background mortality rate

%% Ambient temperature (Kelvin)

T = 273.15 + 20;

%% Temperature- or body mass-dependent parameters
%  Without size-temperature interaction

% % Resource supply density
% RSmax    = 0.0042 * exp( 0.151/(0.00008617*T) );
% RLmax    = RSmax;
% 
% % Consumer maximum ingestion rate
% ICSRSmax = (19 * (B_CS^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_CS;
% ICSRLmax = (19 * (B_CS^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_CS;
% ICLRSmax = (19 * (B_CL^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_CL;
% ICLRLmax = (19 * (B_CL^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_CL;
% 
% % Consumer metabolic rate
% mCS      = (850000000 * (B_CS^0.7) * exp( -0.56/(0.00008617*T) )) / B_CS;
% mCL      = (850000000 * (B_CL^0.7) * exp( -0.56/(0.00008617*T) )) / B_CL;

%% Temperature- or body mass-dependent parameters
%  With size-temperature interaction in Rmax and in temperature optimum of Imax

% Resource supply density
RSmax    = 0.0042 * exp( 0.151/(0.00008617*T) );
RLmax    = (5.88* 10^(-7)) * exp( 0.37564/(0.00008617*T) );

% Consumer maximum ingestion rate
ICSRSmax = (19 * (B_CS^(0.7)) * exp(-((T-(273.15+24))^2)/(2*(8^2)))) / B_CS;
ICSRLmax = (19 * (B_CS^(0.7)) * exp(-((T-(273.15+24))^2)/(2*(8^2)))) / B_CS;
ICLRSmax = (19 * (B_CL^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_CL;
ICLRLmax = (19 * (B_CL^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_CL;

% Consumer metabolic rate
mCS      = (850000000 * (B_CS^0.7) * exp( -0.56/(0.00008617*T) )) / B_CS;
mCL      = (850000000 * (B_CL^0.7) * exp( -0.56/(0.00008617*T) )) / B_CL;

%% Temperature- or body mass-dependent parameters
%  With size-temperature interaction in Rmax and in metabolic rate

% % Resource supply density
% RSmax    = 0.0042 * exp( 0.151/(0.00008617*T) );
% RLmax    = (5.88* 10^(-7)) * exp( 0.37564/(0.00008617*T) );
% 
% % Consumer maximum ingestion rate
% ICSRSmax = (19 * (B_CS^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_CS;
% ICSRLmax = (19 * (B_CS^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_CS;
% ICLRSmax = (19 * (B_CL^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_CL;
% ICLRLmax = (19 * (B_CL^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_CL;
% 
% % Consumer metabolic rate
% mCS      = (850000000 * (B_CS^(0.7 + 0.0005*T)) * exp( -0.56/(0.00008617*T) )) / B_CS;
% mCL      = (850000000 * (B_CL^(0.7 + 0.0005*T)) * exp( -0.56/(0.00008617*T) )) / B_CL;

%% Variables

RS = Y(1);  % small resource biomass density
RL = Y(2);  % large resource biomass density
CS = Y(3);  % small consumer biomass density
CL = Y(4);  % large consumer biomass density

%% Ingestion rates

ICSRS = ( pCSRS * (ICSRSmax/(HCSRS^(1+q))) * RS^(1+q) + 0 * (ICSRLmax/(HCSRL^(1+q))) * RL^(1+q) ) / ... 
    ( 1 + (pCSRS/(HCSRS^(1+q))) * RS^(1+q) + (pCSRL/(HCSRL^(1+q))) * RL^(1+q) );

ICSRL = ( 0 * (ICSRSmax/(HCSRS^(1+q))) * RS^(1+q) + pCSRL * (ICSRLmax/(HCSRL^(1+q))) * RL^(1+q) ) / ... 
    ( 1 + (pCSRS/(HCSRS^(1+q))) * RS^(1+q) + (pCSRL/(HCSRL^(1+q))) * RL^(1+q) );

ICLRS = ( pCLRS * (ICLRSmax/(HCLRS^(1+q))) * RS^(1+q) + 0 * (ICLRLmax/(HCLRL^(1+q))) * RL^(1+q) ) / ... 
    ( 1 + (pCLRS/(HCLRS^(1+q))) * RS^(1+q) + (pCLRL/(HCLRL^(1+q))) * RL^(1+q) );

ICLRL = ( 0 * (ICLRSmax/(HCLRS^(1+q))) * RS^(1+q) + pCLRL * (ICLRLmax/(HCLRL^(1+q))) * RL^(1+q) ) / ... 
    ( 1 + (pCLRS/(HCLRS^(1+q))) * RS^(1+q) + (pCLRL/(HCLRL^(1+q))) * RL^(1+q) );

%% ODE system

dRSdt = deltaRS*(RSmax - RS) - ICSRS*CS - ICLRS*CL;

dRLdt = deltaRL*(RLmax - RL) - ICSRL*CS - ICLRL*CL;

dCSdt =  betaCS*(ICSRS+ICSRL)*CS - mCS*CS - muCS*CS;

dCLdt =  betaCL*(ICLRS+ICLRL)*CL - mCL*CL - muCL*CL;

dYdt = [dRSdt; dRLdt; dCSdt; dCLdt];

end