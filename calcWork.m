function [work, power, theta] = calcWork(force,pos,dt)
%% CALCWORK: Calculates mechanical work
%
%   work = calcWork(force, position, dt) calculates the mechanical work
%   WORK given a 3D force vector FORCE, a 3D position vector POSITION, and 
%   scalar sampling interval DT. FORCE and POSITION are 3-by-N arrays,
%   where N is the number of time points
%
%   [work, power, theta] = calcWork(__) also returns the power transfer
%   POWER and the angle between the force and velocity THETA.
%
%   CALCWORK calculates work by first calculating the instantaneous power
%   transfer, and then integrating the power transfer to estimate the work.

%   Luke Drnach
%   December 4, 2018

%   Slight changes MW 11/22/19

%Calculate instantaneous power as the product of force and velocity
vel = diff(pos)./dt;               % Estimate the velocity
power = [0,dot(force(2:end,:)',vel')]';
% Integrate power to get cumulative work
work = cumsum(power).*dt;
% Calculate the angle between the force and the velocity
theta = power(2:end)./(sqrt(sum(force(2:end,:).^2,2)).*sqrt(sum(vel.^2,2)));
end