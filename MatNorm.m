function dataOut = MatNorm(dataIn,dim)
% �Գ���ľ���dataIn���й�һ��
% ���룺
%   dataIn -- �����ά����
%   dim -- ��һ��ά��: dim = 'col'�� 1, ��ά��Ĭ�ϣ�
%                      dim = 'row'�� 2, ��ά
% �����
%   dataOut -- ��һ����ľ���
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