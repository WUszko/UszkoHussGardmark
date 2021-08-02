%% Community III
% Two stage-structured consumer species feeding on two resources
% For units and references, see Table S1.2 in Appendix S1
% Created by Wojciech Uszko (2021)

function dYdt = model_III(t, Y)

%% Body masses (ng dry weight)

B_JS = 100;     % juvenile of small consumer
B_AS = 1000;    % adult of small consumer
B_JL = 1000;    % juvenile of large consumer
B_AL = 10000;   % adult of large consumer

%% Temperature- or body mass-independent parameters

deltaRS = 0.1;  % small resource supply rate
deltaRL = 0.1;  % large resource supply rate

q       = 0;    % functional response (Hill) exponent; if =0 then type II

p       = 0.75; % diet preference
pASRS   = p;
pASRL   = 1-pASRS;
pJLRL   = pASRS;
pJLRS   = 1-pASRS;

betaJS   = 0.6;  % small juvenile conversion efficiency
betaAS   = 0.6;  % small adult conversion efficiency
betaJL   = 0.6;  % large juvenile conversion efficiency
betaAL   = 0.6;  % large adult conversion efficiency

HJSRS   = 0.2;  % half-saturation constant
HASRS   = 0.2;  % half-saturation constant
HASRL   = 0.2;  % half-saturation constant
HJLRS   = 0.2;  % half-saturation constant
HJLRL   = 0.2;  % half-saturation constant
HALRL   = 0.2;  % half-saturation constant

zJSAS   = 0.1;  % juvenile-to-adult mass ratio in small consumer
zJLAL   = 0.1;  % juvenile-to-adult mass ratio in large consumer

muJS    = 0.01; % small juvenile background mortality rate
muAS    = 0.01; % small adult background mortality rate
muJL    = 0.01; % large juvenile background mortality rate
muAL    = 0.01; % large adult background mortality rate

%% Ambient temperature (Kelvin)

T = 273.15 + 20;

%% Temperature- or body mass-dependent parameters
%  Without size-temperature interaction

% % Resource supply density
% RSmax   = 0.0042 * exp( 0.151/(0.00008617*T) );
% RLmax   = RSmax;
% 
% % Consumer maximum ingestion rate
% IJSRSmax = (19 * (B_JS^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_JS;
% IASRSmax = (19 * (B_AS^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_AS;
% IASRLmax = (19 * (B_AS^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_AS;
% IJLRSmax = (19 * (B_JL^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_JL;
% IJLRLmax = (19 * (B_JL^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_JL;
% IALRLmax = (19 * (B_AL^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_AL;
% 
% % Consumer metabolic rate
% mJS      = (850000000 * (B_JS^0.7) * exp( -0.56/(0.00008617*T) )) / B_JS;
% mAS      = (850000000 * (B_AS^0.7) * exp( -0.56/(0.00008617*T) )) / B_AS;
% mJL      = (850000000 * (B_JL^0.7) * exp( -0.56/(0.00008617*T) )) / B_JL;
% mAL      = (850000000 * (B_AL^0.7) * exp( -0.56/(0.00008617*T) )) / B_AL;

%% Temperature- or body mass-dependent parameters
%  With size-temperature interaction in Rmax and in temperature optimum of Imax

% Resource supply density
RSmax    = 0.0042 * exp( 0.151/(0.00008617*T) );
RLmax    = (5.88* 10^(-7)) * exp( 0.37564/(0.00008617*T) );

% Consumer maximum ingestion rate
IJSRSmax = (19 * (B_JS^(0.7)) * exp(-((T-(273.15+24))^2)/(2*(8^2)))) / B_JS;
IASRSmax = (19 * (B_AS^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_AS;
IASRLmax = (19 * (B_AS^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_AS;
IJLRSmax = (19 * (B_JL^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_JL;
IJLRLmax = (19 * (B_JL^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_JL;
IALRLmax = (19 * (B_AL^(0.7)) * exp(-((T-(273.15+16))^2)/(2*(8^2)))) / B_AL;

% Consumer metabolic rate
mJS      = (850000000 * (B_JS^0.7) * exp( -0.56/(0.00008617*T) )) / B_JS;
mAS      = (850000000 * (B_AS^0.7) * exp( -0.56/(0.00008617*T) )) / B_AS;
mJL      = (850000000 * (B_JL^0.7) * exp( -0.56/(0.00008617*T) )) / B_JL;
mAL      = (850000000 * (B_AL^0.7) * exp( -0.56/(0.00008617*T) )) / B_AL;

%% Temperature- or body mass-dependent parameters
%  With size-temperature interaction in Rmax and in metabolic rate

% % Resource supply density
% RSmax   = 0.0042 * exp( 0.151/(0.00008617*T) );
% RLmax   = (5.88* 10^(-7)) * exp( 0.37564/(0.00008617*T) );
% 
% % Consumer maximum ingestion rate
% IJSRSmax = (19 * (B_JS^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_JS;
% IASRSmax = (19 * (B_AS^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_AS;
% IASRLmax = (19 * (B_AS^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_AS;
% IJLRSmax = (19 * (B_JL^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_JL;
% IJLRLmax = (19 * (B_JL^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_JL;
% IALRLmax = (19 * (B_AL^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_AL;
% 
% % Consumer metabolic rate
% mJS      = (850000000 * (B_JS^(0.7 + 0.0005*T)) * exp( -0.56/(0.00008617*T) )) / B_JS;
% mAS      = (850000000 * (B_AS^(0.7 + 0.0005*T)) * exp( -0.56/(0.00008617*T) )) / B_AS;
% mJL      = (850000000 * (B_JL^(0.7 + 0.0005*T)) * exp( -0.56/(0.00008617*T) )) / B_JL;
% mAL      = (850000000 * (B_AL^(0.7 + 0.0005*T)) * exp( -0.56/(0.00008617*T) )) / B_AL;

%% Variables

RS = Y(1);  % small resource biomass density
RL = Y(2);  % large resource biomass density
JS = Y(3);  % small juvenile biomass density
AS = Y(4);  % small adult biomass density
JL = Y(5);  % large juvenile biomass density 
AL = Y(6);  % large adult biomass density

%% Ingestion rates

IJSRS = ( 1 * (IJSRSmax/(HJSRS^(1+q))) * RS^(1+q) ) / ... 
    ( 1 + (1/(HJSRS^(1+q))) * RS^(1+q) );

IASRS = ( pASRS * (IASRSmax/(HASRS^(1+q))) * RS^(1+q) + 0 * (IASRLmax/(HASRL^(1+q))) * RL^(1+q) ) / ... 
    ( 1 + (pASRS/(HASRS^(1+q))) * RS^(1+q) + (pASRL/(HASRL^(1+q))) * RL^(1+q) );

IASRL = ( 0 * (IASRSmax/(HASRS^(1+q))) * RS^(1+q) + pASRL * (IASRLmax/(HASRL^(1+q))) * RL^(1+q) ) / ... 
    ( 1 + (pASRS/(HASRS^(1+q))) * RS^(1+q) + (pASRL/(HASRL^(1+q))) * RL^(1+q) );

IJLRS = ( pJLRS * (IJLRSmax/(HJLRS^(1+q))) * RS^(1+q) + 0 * (IJLRLmax/(HJLRL^(1+q))) * RL^(1+q) ) / ... 
    ( 1 + (pJLRS/(HJLRS^(1+q))) * RS^(1+q) + (pJLRL/(HJLRL^(1+q))) * RL^(1+q) );

IJLRL = ( 0 * (IJLRSmax/(HJLRS^(1+q))) * RS^(1+q) + pJLRL * (IJLRLmax/(HJLRL^(1+q))) * RL^(1+q) ) / ... 
    ( 1 + (pJLRS/(HJLRS^(1+q))) * RS^(1+q) + (pJLRL/(HJLRL^(1+q))) * RL^(1+q) );

IALRL = ( 1 * (IALRLmax/(HALRL^(1+q))) * RL^(1+q) ) / ... 
    ( 1 + (1/(HALRL^(1+q))) * RL^(1+q) );

%% Stage-specific production rates

vJS     = betaJS*IJSRS - mJS;                           % small juvenile production rate
gammaJS = (vJS - muJS) / (1 - zJSAS^(1-(muJS/vJS)));    % maturation rate
vAS     = betaAS*(IASRS+IASRL) - mAS;                   % reproduction rate

vJL     = betaJL*(IJLRS+IJLRL) - mJL;                   % large juvenile production rate
gammaJL = (vJL - muJL) / (1 - zJLAL^(1-(muJL/vJL)));    % maturation rate
vAL     = betaAL*IALRL - mAL;                           % reproduction rate

%% ODE system

dRSdt = deltaRS*(RSmax - RS) - IJSRS*JS - IASRS*AS - IJLRS*JL;

dRLdt = deltaRL*(RLmax - RL) - IASRL*AS - IJLRL*JL - IALRL*AL;

dJSdt = max(vAS,0)*AS + vJS*JS - max(gammaJS,0)*JS - muJS*JS;

dASdt = max(gammaJS,0)*JS + min(vAS,0)*AS - muAS*AS;

dJLdt = max(vAL,0)*AL + vJL*JL - max(gammaJL,0)*JL - muJL*JL;

dALdt = max(gammaJL,0)*JL + min(vAL,0)*AL - muAL*AL;

dYdt = [dRSdt; dRLdt; dJSdt; dASdt; dJLdt; dALdt];

end