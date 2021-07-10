% ==========================================================================
% Using the 'whos' command and evaluating inside the base workspace and 
%   sumps up the bytes.
% The output is display in MB.
% ==========================================================================
function [ memory_in_use ] = ShowMemory()

    mem_elements = evalin('base','whos');
    if size(mem_elements,1) > 0

        for i = 1:size(mem_elements,1)
            memory_array(i) = mem_elements(i).bytes;
        end

        memory_in_use = sum(memory_array);
        memory_in_use = memory_in_use/1024/1024;
    else
        memory_in_use = 0;
    end

end
