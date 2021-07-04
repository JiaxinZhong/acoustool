% 判断所有输入的变量是否是compatible size array
% 输入：任意数目的数值变量

function result = isCompatibleSize(varargin)

    tmp = 0;
    try 
        for i = 1:length(varargin)
            tmp = tmp .* varargin{i};
        end
    catch ME
        error([ME.message, '\n'])
    end
    result = 1;
end
