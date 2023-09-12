%% FUNCTION - Sizing Based on Payload Mass (One stage)
function [m0,ms,mp,mu,PI] = sizing(dvrq,mPL,sigma,ve)

    % mass ratio mu - initial mass / final mass
    mu = exp(dvrq/(ve));
    % payload ratio pi - payload mass / inital mass
    PI = (1/mu-sigma)/(1-sigma);

    % initial mass m0
    m0 = mPL/PI;
    % structure mass
    ms = m0*sigma*(1-PI);
    % propullent mass
    mp = (mu-1)/mu*m0;
    