function [r,pos] = corr_Max4Matrix(A,dim,flag_simplify,maxNumVect)
% 计算矩阵A指定维度任意两个相异向量的最大相关系数
% 输入：
% A -- 输入观测矩阵
% dim -- 指定维度: 'col' --列（默认）; 'row'--行
% flag_simplify -- 简化计算的标志位 flag_simplify = 's'，随机提取观测值
% maxNumVect -- 随机提取列数
% 输出:
% r -- 最大相关系数(0~1)
% pos -- 所在位置(i,j)
% Jiwen Geng; 2019/7/31
% Update: 
%  2019/7/31    增加简化操作
if nargin ==1
    dim = 'col';
    flag_simplify = [];
    maxNumVect = [];
    index_customed = 1:size(A,2); 
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
            index_customed = randperm(colNum,maxNumVect);   % 否则随机生成抽样的序列
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
r_max = 0;
if strcmpi(dim,'col')
    for jj = 1:Num_list
        vect1_current = A(:,index_customed(jj));   % M*1
        for kk = jj+1:Num_list
            vect2_current = A(:,index_customed(kk));   % M*1
            r_current = dot(vect1_current,vect2_current)/(norm(vect1_current)*norm(vect2_current));
            if abs(r_current)>r_max
                r_max = abs(r_current);
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
            if abs(r_current)>r_max
                r_max = abs(r_current);
                index1 = index_customed(ii);
                index2 = index_customed(kk);
            end
        end
    end
end

r = r_max;
pos = [index1;index2];

end