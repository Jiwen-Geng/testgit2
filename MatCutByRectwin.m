function dataOut = MatCutByRectwin(dataIn,ratio)
% 对输入的二维矩阵dataIn二维矩形窗截断
% 输入：
%   dataIn -- 输入二维矩阵
%   ratio -- 窗截断的比例，默认为0.5
% 输出：
%   dataOut -- 输出二维矩阵（一般为复信号）
% Jiwen Geng   Date: 2019/08/07
if nargin ==1
    ratio = 0.5;
end
[Na,Nr] = size(dataIn);
Nwin_a = max([trans2even(ratio*Na),2]);
Nwin_r = max([trans2even(ratio*Nr),2]);
index_a = 1 + round(Na/2) + (-Nwin_a/2:Nwin_a/2-1);
index_r = 1 + round(Nr/2) + (-Nwin_r/2:Nwin_r/2-1);
data_2df = ftx(fty(dataIn));
data_2dfWin = zeros(Na,Nr);
data_2dfWin(index_a,index_r) = data_2df(index_a,index_r);
dataOut = iftx(ifty(data_2dfWin));

end