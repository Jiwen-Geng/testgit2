function dataOut = MatCutByRectwin(dataIn,ratio)
% ������Ķ�ά����dataIn��ά���δ��ض�
% ���룺
%   dataIn -- �����ά����
%   ratio -- ���ضϵı�����Ĭ��Ϊ0.5
% �����
%   dataOut -- �����ά����һ��Ϊ���źţ�
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