% 功能：声波的信息
classdef SoundWave < handle
	properties
		freq
        note = nan;
		absorp = nan; % 空气吸声系数，单位Np
		wavnum
        sound_speed = 343;
	end

	properties (Dependent)
		angfreq
        wavlen
	end

	methods 
		function obj = SoundWave(freq)
			obj.freq = freq;
            
            obj.CalWavnum;
            obj.CalAbsorpCoef;
		end

		function angfreq = get.angfreq(obj)
			angfreq = 2*pi*obj.freq;
		end

		function CalWavnum(obj)
			obj.wavnum = obj.angfreq / obj.sound_speed;
        end
        

		function wavlen = get.wavlen(obj)
			wavlen = obj.sound_speed./obj.freq;
        end

        function CalAbsorpCoef(obj, varargin)
            obj.absorp = AbsorpCoeff(obj.freq, varargin{:});
            obj.wavnum = obj.wavnum + 1i*obj.absorp;
		end 

	end
		
end
