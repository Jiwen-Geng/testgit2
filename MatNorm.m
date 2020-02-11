function dataOut = MatNorm(dataIn,dim)
% 对出入的矩阵dataIn进行归一化
% 输入：
%   dataIn -- 输入二维矩阵
%   dim -- 归一化维度: dim = 'col'或 1, 列维（默认）
%                      dim = 'row'或 2, 行维
% 输出：
%   dataOut -- 归一化后的矩阵
% Jiwen Geng      Date: 2019/08/07
if nargin == 1
    dim = 'col';
end
[Na,Nr] = size(dataIn);

if strcmpi(dim,'col')||(dim==1)
    dataOut = dataIn./repmat(sqrt(sum(abs(dataIn).^2,1)),Na,1);
elseif strcmpi(dim,'row')||(dim==2)
    dataOut = dataIn./repmat(sqrt(sum(abs(dataIn).^2,2)),1,Nr);
end

end