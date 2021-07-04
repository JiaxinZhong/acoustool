% =========================================================================
% VERSION INFO
%	Last Modified 	--- 2020-06-25
%	Version 		--- 2.4
% -------------------------------------------------------------------------
% FUNCTION  
%	- Calculate the absorption coefficients by the atmospheric parameters
%		based on the standard ISO 9613-1.
% -------------------------------------------------------------------------
% REFERENCES
%	[1] ISO Technical Committees: Noise. Acoustics — Attenuation of sound 
%		during propagation outdoors — Part 1: Calculation of the absorption
%		of sound by the atmosphere: ISO 9613-1:1993[S]. Geneva: 
%		International Organization for Standardization, 1993.
%	[2]	National Physical Laboratory. NPL Acoustics: Calculation of 
%		absorption of sound by the atmosphere[EB/OL]. [2018-08-08].
%		http://resource.npl.co.uk/acoustics/techguides/absorption/.
% -------------------------------------------------------------------------
% DEPENDENCIES
%	- isCompatibleSize.m
% -------------------------------------------------------------------------
% INPUT
%	freq		- Frequency, in Hertz. 
% OPTIONS
%	humidity 	- relative humidity in percentage
%				- dafault: 60
%	temperature	- Ambient atmospheric temperature, in Celcius
%				- default: 25
%	pressure 	- the atmospheric pressure, in kilopscals
%				- default: 101.325
% NOTES
%	- the dimension of all inputs must be compatible
% -------------------------------------------------------------------------
% OUTPUT
%	alpha_Np	- pure-tone sound attenuation coeff. in Neper per meter, for 
%					atmospheric absorption
%	alpha_dB	- pure-tone sound attenutaion oeffi. in dB per meter, for 
%					atmospheric absorption
% =========================================================================

function [alpha_Np, alpha_dB] = AbsorpCoeff(freq, varargin)
    
	p = inputParser();
	addParameter(p, 'temperature', 25);
	addParameter(p, 'pressure', 101.325);
	addParameter(p, 'humidity', 60);
	parse(p, varargin{:});
	ip = p.Results;

	if ~isCompatibleSize(freq, ip.temperature, ip.pressure, ip.humidity)
		error('All inputs must be compatible in sizes!\n')
	end
    
	% reference air temperature in Kelvins, i.e. 20 Celcius degree
	T0 = 293.15; 
	% the triple-point isotherm temperature (i.e. +0.01 Celsius degree)
	T01 = 273.16; 
	% Ambient atmospheric temperature, in Kelvins
	T = ip.temperature + 273.15; 

	pr = 101.325; % reference ambient atmospheric pressure, in kilopascals
    
    C = -6.8346*(T01./T).^1.261 + 4.6151;
    psat = pr .* 10.^C; % saturation vapour pressure
	% the molar concentration of water vapoupr as a percentage
    h = ip.humidity.*(psat./pr).*(ip.pressure./pr); 
    
    % the oxygen relaxation frequency
    f_rO = ip.pressure./pr .* (24 + 4.04*10^4 * h .* (0.02+h) ./ (0.391+h) );
    % the nitrogen relaxation frequency
    f_rN = ip.pressure./pr .* (T./T0).^(-1/2) .* (9 + 280 * h ...
        .* exp(-4.17.*((T./T0).^(-1/3)-1)));
    
    alpha_Np = freq.^2 .* (1.84*10^(-11) .* pr./ip.pressure .* (T./T0).^(1/2) ...
        + (T./T0)^(-5/2) * (0.01275*exp(-2239.1./T)./(f_rO+freq.^2./f_rO) ...
        + 0.1068*exp(-3352.0./T)./(f_rN+freq.^2./f_rN)) );
	alpha_dB = 20/log(10) * alpha_Np;
end

