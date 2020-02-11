function [r,pos] = corr_Min4Matrix(A,dim,flag_simplify,maxNumVect)
% �������Aָ��ά����������������������С���ϵ��
% ���룺
% A -- ����۲����
% dim -- ָ��ά��: 'col' --�У�Ĭ�ϣ�; 'row'--��
% flag_simplify -- �򻯼���ı�־λ flag_simplify = 's'�������ȡ�۲�ֵ
% maxNumVect -- �����ȡ����
% ���:
% r -- ������ϵ��(0~1)
% pos -- ����λ��(i,j)
% Jiwen Geng; 2019/7/31
% Update: 
%  2019/7/31    ���Ӽ򻯲���
if nargin ==1
    dim = 'col';
    flag_simplify = [];
    maxNumVect = [];
elseif nargin ==2
    flag_simplify = [];
    maxNumVect = [];
end
[rowNum,colNum] = size(A);
if strcmpi(flag_simplify,'s')||strcmpi(flag_simplify,'simplify')||...
        strcmpi(flag_simplify,'simplified')||strcmpi(flag_simplify,'less')
    if strcmpi(dim,'col')
        if maxNumVect>= colNum
            index_customed = 1:colNum;
        else
            index_customed = randperm(colNum,maxNumVect);   % ����������ɳ���������
        end
    elseif strcmpi(dim,'row')
        if maxNumVect>= rowNum
            index_customed = 1:rowNum;
        else
            index_customed = randperm(rowNum,maxNumVect);   
        end
    end
end
Num_list = length(index_customed);
r_min = 1;
if strcmpi(dim,'col')
    for jj = 1:Num_list
        vect1_current = A(:,index_customed(jj));   % M*1
        for kk = jj+1:Num_list
            vect2_current = A(:,index_customed(kk));   % M*1
            r_current = dot(vect1_current,vect2_current)/(norm(vect1_current)*norm(vect2_current));
            if abs(r_current)<r_min
                r_min = abs(r_current);
                index1 = index_customed(jj);
                index2 = index_customed(kk);
            end   
        end
    end
elseif strcmpi(dim,'row')
    for ii = 1:Num_list
        vect1_current = A(index_customed(ii),:).';   % N*1
        for kk = ii+1:Num_list
            vect2_current = A(index_customed(kk),:).';   % N*1
            r_current = dot(vect1_current,vect2_current)/(norm(vect1_current)*norm(vect2_current));
            if abs(r_current)<r_min
                r_min = abs(r_current);
                index1 = index_customed(ii);
                index2 = index_customed(kk);
            end
        end
    end
end

r = r_min;
pos = [index1;index2];

end