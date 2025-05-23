function [Q,R,P] = IP_CBD(H,blockNum)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
[m,n] = size(H);
% Q = eye(m);
% P = eye(n);
[U,S,V] = svd(H);
d = diag (S) ;
% d1 = diag (S) ;
l = min (m, n) ;


for len = l : -1 : 1
    if ( d (len) >= eps ) %省略掉小于eps的奇异值 不影响矩阵的数值
        break ;
    end
end

% p = 1:len;
J = len/blockNum;
% pos = zeros(blockNum,T );
% p1 = reshape(p,[blockNum,T]);
% if mod(T,2)==0 %假设偶数根天线 case1 Nt=8或4,二分
%     p1tmp = flip(p1(:,T/2+1:end));
%     pos = [p1(:,1:T/2),p1tmp];
%
% else
%     p1tmp = flip(p1(:,(T+3)/2:end));
%     pos = [p1(:,1:(T+1)/2),p1tmp];
% end
% % pos = reshape(fliplr(pos)',[],1);
% pos = [16 15 2 1  14 13 4 3 12 11 6 5 10 9 8 7];

if mod(J,2)==0 %
    % ps = 1:len/2;
    pss = flip(reshape(1:len/2, J/2,[]),1);
    % pc = len:-1:len/2+1;
    pcs = reshape(len:-1:len/2+1, J/2,[]);
    pos = [pcs.',pss.'];
else
    ps = 1:(len-blockNum)/2;
    pss = flip(reshape(ps, (J-1)/2,[]),1);
    pc = len:-1:(len+blockNum)/2+1;
    pcs = reshape(pc, (J-1)/2,[]);
    pp = (len+blockNum)/2:-1:(len-blockNum)/2+1;
    pos = [pcs.',pp.',pss.'];
end
pos = reshape(pos',[],1);
if blockNum==1
    pos = 1:len;
end
E = eye(n);
E = E(:,pos);
U = U*E;
V = V*E;
S = E'*S*E;
for bb = 1:blockNum
    for i = len/blockNum*(bb-1)+1:len/blockNum*bb-1

        d = diag(S);
        perm = [i i+1];
        delta1 = d (i) ;
        delta2 = d (i+1) ;
        if delta1==delta2
            continue
        end
        t = delta1 + delta2 ;

        sigma_bar = prod(d(perm))^(1/2);
        f = (delta1 - sigma_bar)/(delta1 - delta2) ;
        s = sqrt (f*(delta1+sigma_bar)/t) ;
        c = sqrt ((1-f)*(delta2+sigma_bar)/t) ;

        % d (i+1) = delta1*delta2/sigma_bar ;% = y in paper
        % z (i) = s*c*(delta2 - delta1)*t/sigma_bar ;% = x in paper
        G1 = [ c -s
            s  c ] ;
        V(:,perm) = V(:,perm)*G1;
        G2 = (1/sigma_bar)*[ c*delta1 -s*delta2
            s*delta2  c*delta1 ] ;
        U(:,perm) = U(:,perm)*G2;
        S(perm,:) = G2'*S(perm,:);
        S(:,perm) = S(:,perm)*G1;



        for j=i:-1:len/blockNum*(bb-1)+2
            perm = [j-1 j];
            if S(j-1,j+1)==0 && S(j,j+1)==0
                continue
            end
            Gr = givensB(S(j-1,j+1),S(j,j+1),'ColGivens-d');
            S(perm,:) = Gr*S(perm,:);
            U(:,perm) = U(:,perm)*Gr';
            if S(j,j-1)~=0
                % perm = [i+2 i+1];
                Gl = givensB(S(j,j-1),S(j,j),'RowGivens-d');
                S(:,perm) = S(:,perm)*Gl;
                V(:,perm) = V(:,perm)*Gl;
            end
        end
    end   % 从前向中间迭代
end


Q = U;
P = V;
R = S;
end