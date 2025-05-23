function LLRext = MIMO_block_optimal_detection(yk, Leck, noiseVar2, bidiagMat, M, Nr, BlockNum)
%% ===========================================
% Function: Equalization forward/backward algorithm
% =================================================
% this algorithm is designed for sequence detection
% Input :
%   yk：receive signal；
%   Leck<=>prior probability；
%   noiseVar2 ： noise variance
%   bidiagMat ： bidiagonal matrix
%   M ： modulation order
%   Nt： sender antenna
%   BlockNum： PL-CBD BLOCK num
% Output :
%   LLRext：log lokelihhos ratio LLR
%% ===================================
Qm = log2(M);
noiseVar = noiseVar2;
J = Nr/BlockNum;
if M ==2 %BPSK
    [m,n] = size(yk);
    Lcky = zeros(m,n);
    % state = qammod([1;0], M, 'InputType', 'bit', 'UnitAveragePower', true);
    state = nrSymbolModulate([1; ...
        0;],'BPSK');
    Apos = [1 0 ;
        1 0 ]; %  A(+1)
    Aneg = [0 1 ;
        0 1 ]; %  A(-1)
    for nn = 1:n
        y = yk(:,nn);
        Leckp = Leck(Nr*Qm*(nn-1)+1:Nr*Qm*nn);
        trellisOutput = zeros(M,M,Nr);
        %% 计算状态转移矩阵
        for ii = Nr:-1:1
            if ii == Nr
                for mm = 1:M
                    trellisOutput(:,mm,ii) = bidiagMat(ii,ii)*state(mm); %第Nt根天线的状态图输出
                end
            else
                for mm = 1:M
                    trellisOutput(:,mm,ii) = bidiagMat(ii,ii+1)*state + bidiagMat(ii,ii)*state(mm); %第Nt根天线的状态图输出
                end
            end
        end
        %% 前后向迭代计算
        bk1 =  ones(M,1,Nr+1);  %bk
        for tt = Nr:-1:2   %这里的Nt理解为状态图从后往前
            bk1(:,:,tt) = exp(Apos.*Leckp(Nr-tt+1))/(1 + exp(Leckp(Nr-tt+1)))...
                .*(exp(-(abs(y(Nr-tt+1)-trellisOutput(:,:,Nr-tt+1))).^2/noiseVar))*bk1(:,:,tt+1);
            bk1(:,:,tt) = bk1(:,:,tt)/sum(bk1(:,:,tt));
        end
        fk = ones(M,1); %f0


        for kk = 1 : Nr

            % gammakprior = exp(Apos.*Leckp(Nt-kk+2))/(1 + exp(Leckp(Nt-kk+2))); % P(xk=xi,j)
            % PkPrior = gammakprior.*(exp(-(abs(y(Nt-kk+2)-trellisOutput(:,:,Nt-kk+2))).^2/noiseVar));  % P,gamma构成的P矩阵

            gammak = exp(Apos.*Leckp(Nr-kk+1))/(1 + exp(Leckp(Nr-kk+1)));
            Pk=gammak.*(exp(-(abs(y(Nr-kk+1)-trellisOutput(:,:,Nr-kk+1))).^2/noiseVar));
            Lcky(Nr-kk+1,nn) = log(fk'*(Apos.*Pk)*bk1(:,:,kk+1)/(fk'*(Aneg.*Pk)*bk1(:,:,kk+1)));
            fk = Pk'*fk;   % 前向更新 fk
            fk = fk / sum(fk);

        end
    end
    LLRext = reshape(Lcky, size(Leck));
elseif M == 4 % QPSK
    %% ================================
    %-1 1 [1;0;]
    %% ================================
    [m,n] = size(yk);
    Lcky = zeros(Qm*m,n);
    % Lckytmp = zeros(Qm*m,n);
        
    RealQm = Qm/2;
    state = (-1:2:1)'/sqrt(2);

    Afirst1 =  [ones(2,1),zeros(2,1)]; %  A(+1,***)
    Afirst0 =  ones(2)-Afirst1; %  A(-1,***)

    %% bit 间隔排列
    % noiseVar = noiseVar2/2;
    %% 计算状态转移矩阵  虚部实部共用一个，因为信道和状态都是给定的
    trellisOutputAll= ones(sqrt(M),sqrt(M),Nr);
    for ii = Nr:-1:1
        if ii == Nr
            for mm = 1:sqrt(M)
                trellisOutputAll(:,mm,ii) = bidiagMat(ii,ii)*state(mm); %第Nt根天线的状态图输出
            end
        else
            for mm = 1:sqrt(M)
                trellisOutputAll(:,mm,ii) = bidiagMat(ii,ii+1)*state + bidiagMat(ii,ii)*state(mm); %第Nt根天线的状态图输出
            end
        end
    end
    for nn = 1:n

        % Leckp = Leck(Qm*Nt*(nn-1)+1:Qm*Nt*nn);
        Realy = real(yk(:,nn));
        Imagy = imag(yk(:,nn));

        LeckpAll = Leck(Nr*Qm*(nn-1)+1:Nr*Qm*nn);
        Lckytmp = zeros(2,RealQm*m); % LLR of real store in row:1;LLR of imag store in row:2

        %% ================================== start
        for bn = 1:BlockNum
            Leckp = LeckpAll(J*Qm*(bn-1)+1:J*Qm*bn);
            LckytmpBN = zeros(2,RealQm*J);
            RealyBN = Realy(J*(bn-1)+1:J*bn);
            ImagyBN = Imagy(J*(bn-1)+1:J*bn);
            trellisOutput = trellisOutputAll(:,:,J*(bn-1)+1:J*bn);
            %% real
            % 后向传递状态信息
            LeReal = Leckp(1:2:end);

            bkReal =  ones(sqrt(M),1,J+1);  %b_Nt+1
            for tt = J:-1:2   %这里的Nt理解为状态图从后往前
                gammakpost = exp(Afirst1.*LeReal(Qm/2*(J-tt+1)-(Qm/2-1)))/(1 + exp(LeReal(Qm/2*(J-tt+1)-(Qm/2-1))));
                PkPost = gammakpost.* (exp(-(abs(RealyBN(J-tt+1)-trellisOutput(:,:,J-tt+1))).^2/noiseVar)); % Pk,gamma构成的P矩阵

                bkReal(:,:,tt) = PkPost*bkReal(:,:,tt+1);
                bkReal(:,:,tt) = bkReal(:,:,tt)/sum(bkReal(:,:,tt));
            end
            fkReal = ones(sqrt(M),1,J); %f0
            for kk = 1 : J

                gammak = exp(Afirst1.*LeReal(Qm/2*(J-kk+1)-(Qm/2-1)))/(1 + exp(LeReal(Qm/2*(J-kk+1)-(Qm/2-1))));
                Pk=gammak.*exp(-(abs(RealyBN(J-kk+1)-trellisOutput(:,:,J-kk+1))).^2/noiseVar);
                LckytmpBN(1,RealQm*(J-kk)+1) = log(fkReal(:,:,kk)'*(Afirst1.*Pk)*bkReal(:,:,kk+1)/(fkReal(:,:,kk)'*(Afirst0.*Pk)*bkReal(:,:,kk+1)));

                fkReal(:,:,kk+1) = Pk'*fkReal(:,:,kk);   % 前向更新 fk
                fkReal(:,:,kk+1) = fkReal(:,:,kk+1) / sum(fkReal(:,:,kk+1));
            end

            %% imag
            Leimag = Leckp(2:2:end);
            bkImag =  ones(sqrt(M),1,J+1);  %b_Nt+1
            for tt = J:-1:2   %这里的Nt理解为状态图从后往前
                gammakpost = exp(Afirst1.*Leimag(Qm/2*(J-tt+1)-(Qm/2-1)))/(1 + exp(Leimag(Qm/2*(J-tt+1)-(Qm/2-1))));
                PkPost = gammakpost.* (exp(-(abs(ImagyBN(J-tt+1)-trellisOutput(:,:,J-tt+1))).^2/noiseVar)); % Pk,gamma构成的P矩阵

                bkImag(:,:,tt) = PkPost*bkImag(:,:,tt+1);
                bkImag(:,:,tt) = bkImag(:,:,tt)/sum(bkImag(:,:,tt));
            end
            fkImag = ones(sqrt(M),1,J); %f0
            for kk = 1 : J

                gammak = exp(Afirst1.*Leimag(Qm/2*(J-kk+1)-(Qm/2-1)))/(1 + exp(Leimag(Qm/2*(J-kk+1)-(Qm/2-1))));

                Pk=gammak.*exp(-(abs(ImagyBN(J-kk+1)-trellisOutput(:,:,J-kk+1))).^2/noiseVar);
                LckytmpBN(2,RealQm*(J-kk)+1) = log(fkImag(:,:,kk)'*(Afirst1.*Pk)*bkImag(:,:,kk+1)/(fkImag(:,:,kk)'*(Afirst0.*Pk)*bkImag(:,:,kk+1)));

                fkImag(:,:,kk+1) = Pk'*fkImag(:,:,kk);   % 前向更新 fk
                fkImag(:,:,kk+1) = fkImag(:,:,kk+1) / sum(fkImag(:,:,kk+1));
            end
            %%  test end
            Lckytmp(:,RealQm*(J*(bn-1)+1:J*bn)) = LckytmpBN;
        end
        Lcky(:,nn) = reshape(Lckytmp,[],1);
        %% ================================== end
        % %% real
        % % 后向传递状态信息
        % bkReal =  ones(sqrt(M),1,Nt+1);  %b_Nt+1
        % LeReal = Leckp(1:2:end);
        % for tt = Nt:-1:2   %这里的Nt理解为状态图从后往前
        %     gammakpost = exp(Afirst1.*LeReal(RealQm*(Nt-tt+1)-(RealQm-1)))/(1 + exp(LeReal(RealQm*(Nt-tt+1)-(RealQm-1))));
        %     PkPost = gammakpost.* (exp(-((Realy(Nt-tt+1)-trellisOutput(:,:,Nt-tt+1))).^2/noiseVar)); % Pk,gamma构成的P矩阵
        %
        %     bkReal(:,:,tt) = PkPost*bkReal(:,:,tt+1);
        %     bkReal(:,:,tt) = bkReal(:,:,tt)/sum(bkReal(:,:,tt));
        % end
        % fkReal = ones(sqrt(M),1,Nt); %f0
        % for kk = 1 : Nt
        %
        %     gammak = exp(Afirst1.*LeReal(RealQm*(Nt-kk+1)-(RealQm-1)))/(1 + exp(LeReal(RealQm*(Nt-kk+1)-(RealQm-1))));
        %
        %     Pk=gammak.*exp(-((Realy(Nt-kk+1)-trellisOutput(:,:,Nt-kk+1))).^2/noiseVar);
        %     Lckytmp(1,RealQm*(Nt-kk)+1) = log(fkReal(:,:,kk)'*(Afirst1.*Pk)*bkReal(:,:,kk+1)/(fkReal(:,:,kk)'*(Afirst0.*Pk)*bkReal(:,:,kk+1)));
        %
        %     fkReal(:,:,kk+1) = Pk'*fkReal(:,:,kk);   % 前向更新 fk
        %     fkReal(:,:,kk+1) = fkReal(:,:,kk+1) / sum(fkReal(:,:,kk+1));
        % end
        %
        % %% imag
        % Leimag = Leckp(2:2:end);
        % bkImag =  ones(sqrt(M),1,Nt+1);  %b_Nt+1
        % for tt = Nt:-1:2   %这里的Nt理解为状态图从后往前
        %     gammakpost = exp(Afirst1.*Leimag(RealQm*(Nt-tt+1)-(RealQm-1)))/(1 + exp(Leimag(RealQm*(Nt-tt+1)-(RealQm-1))));
        %     PkPost = gammakpost.* (exp(-((Imagy(Nt-tt+1)-trellisOutput(:,:,Nt-tt+1))).^2/noiseVar)); % Pk,gamma构成的P矩阵
        %
        %     bkImag(:,:,tt) = PkPost*bkImag(:,:,tt+1);
        %     bkImag(:,:,tt) = bkImag(:,:,tt)/sum(bkImag(:,:,tt));
        % end
        % fkImag = ones(sqrt(M),1,Nt); %f0
        % for kk = 1 : Nt
        %
        %     gammak = exp(Afirst1.*Leimag(RealQm*(Nt-kk+1)-(RealQm-1)))/(1 + exp(Leimag(RealQm*(Nt-kk+1)-(RealQm-1))));
        %
        %     Pk=gammak.*exp(-((Imagy(Nt-kk+1)-trellisOutput(:,:,Nt-kk+1))).^2/noiseVar);
        %     Lckytmp(2,RealQm*(Nt-kk)+1) = log(fkImag(:,:,kk)'*(Afirst1.*Pk)*bkImag(:,:,kk+1)/(fkImag(:,:,kk)'*(Afirst0.*Pk)*bkImag(:,:,kk+1)));
        %
        %     fkImag(:,:,kk+1) = Pk'*fkImag(:,:,kk);   % 前向更新 fk
        %     fkImag(:,:,kk+1) = fkImag(:,:,kk+1) / sum(fkImag(:,:,kk+1));
        % end
        % %%  test end
        % Lcky(:,nn) = reshape(Lckytmp,[],1);
    end
    LLRext = reshape(Lcky, size(Leck));
elseif M == 16 % 16QAM
    %% ================================
    %-3 -1 [1;1; 1;0;]
    % 1  3 [0;0; 0;1;]
    %% ================================
    [m,n] = size(yk);
    Lcky = zeros(Qm*m,n);
    % Lckytmp = zeros(Qm*m,n);

    RealQm = Qm/2;
    state = (-3:2:3)'/sqrt(10);

    Afirst1 =  [ones(4,2),zeros(4,2)]; %  A(+1,***)
    Afirst0 =  ones(4)-Afirst1; %  A(-1,***)

    Asecond1 =  [ones(4,1),zeros(4,2),ones(4,1)]; %  A(*+1**)
    Asecond0 =  ones(4)-Asecond1;        %  A(*-1,**)

    trellisOutputAll= ones(sqrt(M),sqrt(M),Nr);
    for ii = Nr:-1:1
        if ii == Nr
            for mm = 1:sqrt(M)
                trellisOutputAll(:,mm,ii) = bidiagMat(ii,ii)*state(mm); %第Nt根天线的状态图输出
            end
        else
            for mm = 1:sqrt(M)
                trellisOutputAll(:,mm,ii) = bidiagMat(ii,ii+1)*state + bidiagMat(ii,ii)*state(mm); %第Nt根天线的状态图输出
            end
        end
    end
    %% bit 间隔排列

    for nn = 1:n
        % noiseVar = noiseVar2/2;
        LeckpAll = Leck(Qm*Nr*(nn-1)+1:Qm*Nr*nn);
        Realy = real(yk(:,nn));
        Imagy = imag(yk(:,nn));
        % trellisOutput= ones(sqrt(M),sqrt(M),Nt);
        Lckytmp = zeros(2,RealQm*m);

        %% ================================== start
        for bn = 1:BlockNum
            Leckp = LeckpAll(J*Qm*(bn-1)+1:J*Qm*bn);
            LckytmpBN = zeros(2,RealQm*J);
            RealyBN = Realy(J*(bn-1)+1:J*bn);
            ImagyBN = Imagy(J*(bn-1)+1:J*bn);
            trellisOutput = trellisOutputAll(:,:,J*(bn-1)+1:J*bn);
            %% real
            % 后向传递状态信息
            LeReal = Leckp(1:2:end);
            %% real
            % 后向传递状态信息
            bkReal =  ones(sqrt(M),1,J+1);  %b_Nt+1
            for tt = J:-1:2   %这里的Nt理解为状态图从后往前
                % gammakpost = 1/sqrt(M)*ones(sqrt(M));
                gammakpost =  exp(Afirst1.*LeReal(Qm/2*(J-tt+1)-(Qm/2-1)))/(1 + exp(LeReal(Qm/2*(J-tt+1)-(Qm/2-1)))).*...
                (exp(Asecond1.*LeReal(Qm/2*(J-tt+1)-(Qm/2-2)))/(1 + exp(LeReal(Qm/2*(J-tt+1)-(Qm/2-2)))));
                PkPost = gammakpost.* (exp(-(abs(RealyBN(J-tt+1)-trellisOutput(:,:,J-tt+1))).^2/noiseVar)); % Pk,gamma构成的P矩阵
                
                bkReal(:,:,tt) = PkPost*bkReal(:,:,tt+1);
                if sum(bkReal(:,:,tt))==0
                    bkReal(:,:,tt) = ones(sqrt(M),1);
                end
                bkReal(:,:,tt) = bkReal(:,:,tt)/sum(bkReal(:,:,tt));
            end

            fkReal = ones(sqrt(M),1,J); %f0
            for kk = 1 : J

                % gammak = 1/sqrt(M)*ones(sqrt(M));
                gammak = exp(Afirst1.*LeReal(Qm/2*(J-kk+1)-(Qm/2-1)))/(1 + exp(LeReal(Qm/2*(J-kk+1)-(Qm/2-1)))).*...
                exp(Asecond1.*LeReal(Qm/2*(J-kk+1)-(Qm/2-2)))/(1 + exp(LeReal(Qm/2*(J-kk+1)-(Qm/2-2))));
                Pk=gammak.*exp(-(abs(RealyBN(J-kk+1)-trellisOutput(:,:,J-kk+1))).^2/noiseVar);
                Pk = Pk./sum(sum(Pk));
                X = ((RealyBN(J-kk+1)-trellisOutput(:,:,J-kk+1))).^2;
                [~,c]=find(X==min(min(X)));
                
                LckytmpBN(1,RealQm*(J-kk)+1) = log(fkReal(:,:,kk)'*(Afirst1.*Pk)*bkReal(:,:,kk+1)/(fkReal(:,:,kk)'*(Afirst0.*Pk)*bkReal(:,:,kk+1)));
                if isnan(LckytmpBN(1,RealQm*(J-kk)+1)) && ismember(c(1),[1,2,3,4])
                    LckytmpBN(1,RealQm*(J-kk)+1) = Inf;
                elseif isnan(LckytmpBN(1,RealQm*(J-kk)+1)) && ismember(c(1),[5,6,7,8])
                    LckytmpBN(1,RealQm*(J-kk)+1) = -Inf;
                end
                
                LckytmpBN(1,RealQm*(J-kk)+2) = log(fkReal(:,:,kk)'*(Asecond1.*Pk)*bkReal(:,:,kk+1)/(fkReal(:,:,kk)'*(Asecond0.*Pk)*bkReal(:,:,kk+1)));

                fkReal(:,:,kk+1) = Pk'*fkReal(:,:,kk);   % 前向更新 fk
                fkReal(:,:,kk+1) = fkReal(:,:,kk+1) / sum(fkReal(:,:,kk+1));
            end

            %% imag
            Leimag = Leckp(2:2:end);
            bkImag =  ones(sqrt(M),1,J+1);  %b_Nt+1
            for tt = J:-1:2   %这里的Nt理解为状态图从后往前
                % gammakpost = 1/sqrt(M)*ones(sqrt(M));
                gammakpost = exp(Afirst1.*Leimag(Qm/2*(J-tt+1)-(Qm/2-1)))/(1 + exp(Leimag(Qm/2*(J-tt+1)-(Qm/2-1)))).*...
                (exp(Asecond1.*Leimag(Qm/2*(J-tt+1)-(Qm/2-2)))/(1 + exp(Leimag(Qm/2*(J-tt+1)-(Qm/2-2)))));
                PkPost = gammakpost.* (exp(-(abs(ImagyBN(J-tt+1)-trellisOutput(:,:,J-tt+1))).^2/noiseVar)); % Pk,gamma构成的P矩阵

                bkImag(:,:,tt) = PkPost*bkImag(:,:,tt+1);
                bkImag(:,:,tt) = bkImag(:,:,tt)/sum(bkImag(:,:,tt));
            end
            fkImag = ones(sqrt(M),1,J); %f0
            for kk = 1 : J

                % gammak = 1/sqrt(M)*ones(sqrt(M));
                gammak =exp(Afirst1.*Leimag(Qm/2*(J-kk+1)-(Qm/2-1)))/(1 + exp(Leimag(Qm/2*(J-kk+1)-(Qm/2-1)))).*...
                exp(Asecond1.*Leimag(Qm/2*(J-kk+1)-(Qm/2-2)))/(1 + exp(Leimag(Qm/2*(J-kk+1)-(Qm/2-2))));
                
                Pk=gammak.*exp(-(abs(ImagyBN(J-kk+1)-trellisOutput(:,:,J-kk+1))).^2/noiseVar);
                LckytmpBN(2,RealQm*(J-kk)+1) = log(fkImag(:,:,kk)'*(Afirst1.*Pk)*bkImag(:,:,kk+1)/(fkImag(:,:,kk)'*(Afirst0.*Pk)*bkImag(:,:,kk+1)));
                LckytmpBN(2,RealQm*(J-kk)+2) = log(fkImag(:,:,kk)'*(Asecond1.*Pk)*bkImag(:,:,kk+1)/(fkImag(:,:,kk)'*(Asecond0.*Pk)*bkImag(:,:,kk+1)));

                fkImag(:,:,kk+1) = Pk'*fkImag(:,:,kk);   % 前向更新 fk
                fkImag(:,:,kk+1) = fkImag(:,:,kk+1) / sum(fkImag(:,:,kk+1));
            end
            %%  test end
            Lckytmp(:,RealQm*J*(bn-1)+1:RealQm*J*bn) = LckytmpBN;
        end
        Lcky(:,nn) = reshape(Lckytmp,[],1);
        %% ================================== end
    end
    LLRext = reshape(Lcky, size(Leck));
elseif M== 64
    %% ================================
    %-7:2:-1 [1;1;1; 1;1;0; 1;0;0; 1;0;1;]
    % 1:2:7  [0;0;1; 0;0;0; 0;1;0; 0;1;1;]
    %% ================================
    [m,n] = size(yk);
    Lcky = zeros(Qm*m,n);
    % noiseVar = noiseVar2/2;
    RealQm = Qm/2;
    state = (-7:2:7)'/sqrt(42);
    Afirst1 =  [ones(8,4),zeros(8,4)];
    Afirst0 =  ones(8)-Afirst1;
    Asecond1 =  [ones(8,2),zeros(8,4),ones(8,2)]; %  A(*+1**)
    Asecond0 =  ones(8)-Asecond1;        %  A(*-1,**)
    Athird1 = [ones(8,1),zeros(8,2),ones(8,2),zeros(8,2),ones(8,1)];%  A(**+1*)
    Athird0 = ones(8)-Athird1; %  A(**-1*)

    trellisOutputAll= ones(sqrt(M),sqrt(M),Nr);
    for ii = Nr:-1:1
        if ii == Nr
            for mm = 1:sqrt(M)
                trellisOutputAll(:,mm,ii) = bidiagMat(ii,ii)*state(mm); %第Nt根天线的状态图输出
            end
        else
            for mm = 1:sqrt(M)
                trellisOutputAll(:,mm,ii) = bidiagMat(ii,ii+1)*state + bidiagMat(ii,ii)*state(mm); %第Nt根天线的状态图输出
            end
        end
    end

    for nn = 1:n
        % Leckp = Leck(Qm*Nt*(nn-1)+1:Qm*Nt*nn);
        Realy = real(yk(:,nn));
        Imagy = imag(yk(:,nn));
        Lckytmp = zeros(2,RealQm*m);
        %% ================================== start
        for bn = 1:BlockNum
            LckytmpBN = zeros(2,RealQm*J);
            RealyBN = Realy(J*(bn-1)+1:J*bn);
            ImagyBN = Imagy(J*(bn-1)+1:J*bn);
            trellisOutput = trellisOutputAll(:,:,J*(bn-1)+1:J*bn);
            %% real
            % 后向传递状态信息
            % LeReal = Leckp(1:2:end);
            %% real
            % 后向传递状态信息
            bkReal =  ones(sqrt(M),1,J+1);  %b_Nt+1
            for tt = J:-1:2   %这里的Nt理解为状态图从后往前
                gammakpost = 1/sqrt(M)*ones(sqrt(M));
                PkPost = gammakpost.* (exp(-(abs(RealyBN(J-tt+1)-trellisOutput(:,:,J-tt+1))).^2/noiseVar)); % Pk,gamma构成的P矩阵

                bkReal(:,:,tt) = PkPost*bkReal(:,:,tt+1);
                if sum(bkReal(:,:,tt))==0
                    bkReal(:,:,tt) = ones(sqrt(M),1);
                end
                bkReal(:,:,tt) = bkReal(:,:,tt)/sum(bkReal(:,:,tt));
            end

            fkReal = ones(sqrt(M),1,J); %f0
            for kk = 1 : J

                gammak = 1/sqrt(M)*ones(sqrt(M));
                Pk=gammak.*exp(-(abs(RealyBN(J-kk+1)-trellisOutput(:,:,J-kk+1))).^2/noiseVar);
                Pk = Pk./sum(sum(Pk));
                X = ((RealyBN(J-kk+1)-trellisOutput(:,:,J-kk+1))).^2;
                [~,c]=find(X==min(min(X)));
                
                LckytmpBN(1,RealQm*(J-kk)+1) = log(fkReal(:,:,kk)'*(Afirst1.*Pk)*bkReal(:,:,kk+1)/(fkReal(:,:,kk)'*(Afirst0.*Pk)*bkReal(:,:,kk+1)));
                if isnan(LckytmpBN(1,RealQm*(J-kk)+1)) && ismember(c(1),[1,2,3,4])
                    LckytmpBN(1,RealQm*(J-kk)+1) = Inf;
                elseif isnan(LckytmpBN(1,RealQm*(J-kk)+1)) && ismember(c(1),[5,6,7,8])
                    LckytmpBN(1,RealQm*(J-kk)+1) = -Inf;
                end
                LckytmpBN(1,RealQm*(J-kk)+2) = log(fkReal(:,:,kk)'*(Asecond1.*Pk)*bkReal(:,:,kk+1)/(fkReal(:,:,kk)'*(Asecond0.*Pk)*bkReal(:,:,kk+1)));
                if isnan(LckytmpBN(1,RealQm*(J-kk)+2)) && ismember(c(1),[1,2,7,8])
                    LckytmpBN(1,RealQm*(J-kk)+2) = Inf;
                elseif isnan(LckytmpBN(1,RealQm*(J-kk)+2)) && ismember(c(1),[3,4,5,6])
                    LckytmpBN(1,RealQm*(J-kk)+2) = -Inf;
                end
                LckytmpBN(1,RealQm*(J-kk)+3) = log(fkReal(:,:,kk)'*(Athird1.*Pk)*bkReal(:,:,kk+1)/(fkReal(:,:,kk)'*(Athird0.*Pk)*bkReal(:,:,kk+1)));
                if isnan(LckytmpBN(1,RealQm*(J-kk)+3)) && ismember(c(1),[1,4,5,8])
                    LckytmpBN(1,RealQm*(J-kk)+3) = Inf;
                elseif isnan(LckytmpBN(1,RealQm*(J-kk)+3)) && ismember(c(1),[6,2,3,7])
                    LckytmpBN(1,RealQm*(J-kk)+3) = -Inf;
                end
                fkReal(:,:,kk+1) = Pk'*fkReal(:,:,kk);   % 前向更新 fk
                if sum(fkReal(:,:,kk+1))==0
                    fkReal(:,:,kk+1) = ones(sqrt(M),1);
                end
                fkReal(:,:,kk+1) = fkReal(:,:,kk+1) / sum(fkReal(:,:,kk+1));
            end

            %% imag
            % Leimag = Leckp(2:2:end);
            bkImag =  ones(sqrt(M),1,J+1);  %b_Nt+1
            for tt = J:-1:2   %这里的Nt理解为状态图从后往前
                gammakpost = 1/sqrt(M)*ones(sqrt(M));
                PkPost = gammakpost.* (exp(-(abs(ImagyBN(J-tt+1)-trellisOutput(:,:,J-tt+1))).^2/noiseVar)); % Pk,gamma构成的P矩阵

                bkImag(:,:,tt) = PkPost*bkImag(:,:,tt+1);
                if sum(bkImag(:,:,tt))==0
                    bkImag(:,:,tt) = ones(sqrt(M),1);
                end
                bkImag(:,:,tt) = bkImag(:,:,tt)/sum(bkImag(:,:,tt));
            end
            fkImag = ones(sqrt(M),1,J); %f0
            for kk = 1 : J

                gammak = 1/sqrt(M)*ones(sqrt(M));

                Pk=gammak.*exp(-(abs(ImagyBN(J-kk+1)-trellisOutput(:,:,J-kk+1))).^2/noiseVar);
                Pk = Pk./sum(sum(Pk));
                Y = ((ImagyBN(J-kk+1)-trellisOutput(:,:,J-kk+1))).^2;
                [~,c]=find(Y==min(min(Y)));
                % [~,c]=find(Pk==max(max(Pk)));
                LckytmpBN(2,RealQm*(J-kk)+1) = log(fkImag(:,:,kk)'*(Afirst1.*Pk)*bkImag(:,:,kk+1)/(fkImag(:,:,kk)'*(Afirst0.*Pk)*bkImag(:,:,kk+1)));
                if isnan(LckytmpBN(2,RealQm*(J-kk)+1)) && ismember(c(1),[1,2,3,4])
                    LckytmpBN(2,RealQm*(J-kk)+1) = Inf;
                elseif isnan(LckytmpBN(2,RealQm*(J-kk)+1)) && ismember(c(1),[5,6,7,8])
                    LckytmpBN(2,RealQm*(J-kk)+1) = -Inf;
                end
                LckytmpBN(2,RealQm*(J-kk)+2) = log(fkImag(:,:,kk)'*(Asecond1.*Pk)*bkImag(:,:,kk+1)/(fkImag(:,:,kk)'*(Asecond0.*Pk)*bkImag(:,:,kk+1)));
                if isnan(LckytmpBN(2,RealQm*(J-kk)+2)) && ismember(c(1),[1,2,7,8])
                    LckytmpBN(2,RealQm*(J-kk)+2) = Inf;
                elseif isnan(LckytmpBN(2,RealQm*(J-kk)+2)) && ismember(c(1),[3,4,5,6])
                    LckytmpBN(2,RealQm*(J-kk)+2) = -Inf;
                end
                LckytmpBN(2,RealQm*(J-kk)+3) = log(fkImag(:,:,kk)'*(Athird1.*Pk)*bkImag(:,:,kk+1)/(fkImag(:,:,kk)'*(Athird0.*Pk)*bkImag(:,:,kk+1)));
                if isnan(LckytmpBN(2,RealQm*(J-kk)+3)) && ismember(c(1),[1,4,5,8])
                    LckytmpBN(2,RealQm*(J-kk)+3) = Inf;
                elseif isnan(LckytmpBN(2,RealQm*(J-kk)+3)) && ismember(c(1),[6,2,3,7])
                    LckytmpBN(2,RealQm*(J-kk)+3) = -Inf;
                end
                fkImag(:,:,kk+1) = Pk'*fkImag(:,:,kk);   % 前向更新 fk
                if sum(fkImag(:,:,kk+1))==0
                    fkImag(:,:,kk+1) = ones(sqrt(M),1);
                end
                fkImag(:,:,kk+1) = fkImag(:,:,kk+1) / sum(fkImag(:,:,kk+1));
            end
            %%  test end
            Lckytmp(:,RealQm*J*(bn-1)+1:RealQm*J*bn) = LckytmpBN;
        end
        Lcky(:,nn) = reshape(Lckytmp,[],1);
        %% ================================== end
    end
    LLRext = reshape(Lcky, size(Leck));

elseif M==256
    %% ================================
    %-15:2:-1 [1;1;1;1; 1;1;1;0; 1;1;0;0; 1;1;0;1; 1;0;0;1; 1;0;0;0; 1;0;1;0; 1;0;1;1;]
    % 1:2:15  [0;0;1;1; 0;0;1;0; 0;0;0;0; 0;0;0;1; 0;1;0;1; 0;1;0;0; 0;1;1;0; 0;1;1;1;]
    %% ================================
    [m,n] = size(yk);
    Lcky = zeros(Qm*m,n);
    % Lckytmp = zeros(Qm*m,n);
    % noiseVar = noiseVar2/2;
    RealQm = Qm/2;
    % LckytmpImag = zeros(Qm*m/2,1);

    % Num = 2^Qm;
    % Bit = [0,1];
    % symbolSet = zeros(Qm,M);
    % for nn = 1:Qm
    %     for tt = 1:2^nn
    %         k = mod(tt,2);
    %         if k==0
    %             k=2;
    %         end
    %         symbolSet(nn,2^(Qm-nn)*(tt-1)+1:tt*2^(Qm-nn)) = repmat(Bit(k),[1,2^(Qm-nn)]);
    %     end
    % end  %for循环 给出所有天线发射信号的可能信号矢量symbol
    %
    % state = nrSymbolModulate(reshape(symbolSet,[],1), "256QAM");
    state = (-15:2:15)'/sqrt(170);


    Afirst1 =  [ones(16,8),zeros(16,8)]; %  A(+1,***)
    Afirst0 =  ones(16)-Afirst1; %  A(-1,***)

    Asecond1 =  [ones(16,4),zeros(16,8),ones(16,4)]; %  A(*+1**)
    Asecond0 =  ones(16)-Asecond1;        %  A(*-1,**)

    Athird1 = [ones(16,2),zeros(16,4),ones(16,4),zeros(16,4),ones(16,2)];%  A(**+1*)
    Athird0 = ones(16)-Athird1; %  A(**-1*)

    Aforth1 = [ones(16,1),zeros(16,2),ones(16,2),zeros(16,2),ones(16,2),...
        zeros(16,2),ones(16,2),zeros(16,2),ones(16,1)];%  A(***+1)
    Aforth0 = ones(16)-Aforth1;

    %% bit 间隔排列
    trellisOutputAll= ones(sqrt(M),sqrt(M),Nr);
    for ii = Nr:-1:1  %upper
        if ii == Nr
            for mm = 1:sqrt(M)
                trellisOutputAll(:,mm,ii) = bidiagMat(ii,ii)*state(mm); %第Nt根天线的状态图输出
            end
        else
            for mm = 1:sqrt(M)
                trellisOutputAll(:,mm,ii) = bidiagMat(ii,ii+1)*state + bidiagMat(ii,ii)*state(mm); %第Nt根天线的状态图输出
            end
        end
    end
    %% 计算状态转移矩阵
    % for ii = 1:Nr %lower
    %     if ii == 1
    %         for mm = 1:sqrt(M)
    %             trellisOutputAll(:,mm,ii) = bidiagMat(ii,ii)*state(mm); %第1根天线的状态图输出
    %         end
    %     else
    %         for mm = 1:sqrt(M)
    %             trellisOutputAll(:,mm,ii) = bidiagMat(ii,ii-1)*state + bidiagMat(ii,ii)*state(mm); %第Nt根天线的状态图输出
    %         end
    %     end
    % end
    for nn = 1:n
        % Leckp = Leck(Qm*Nt*(nn-1)+1:Qm*Nt*nn);
        Realy = real(yk(:,nn));
        Imagy = imag(yk(:,nn));
        Lckytmp = zeros(2,RealQm*m);
        

        %% ================================== start
        for bn = 1:BlockNum
            LckytmpBN = zeros(2,RealQm*J);
            RealyBN = Realy(J*(bn-1)+1:J*bn);
            ImagyBN = Imagy(J*(bn-1)+1:J*bn);
            trellisOutput = trellisOutputAll(:,:,J*(bn-1)+1:J*bn);
            %% real
            % 后向传递状态信息
            % LeReal = Leckp(1:2:end);
            %% real
            % 后向传递状态信息
            bkReal =  ones(sqrt(M),1,J+1);  %b_Nt+1
            for tt = J:-1:2   %这里的Nt理解为状态图从后往前
                gammakpost = 1/sqrt(M)*ones(sqrt(M));
                PkPost = gammakpost.* (exp(-((RealyBN(J-tt+1)-trellisOutput(:,:,J-tt+1))).^2/noiseVar)); % Pk,gamma构成的P矩阵

                bkReal(:,:,tt) = PkPost*bkReal(:,:,tt+1);
                bkReal(:,:,tt) = bkReal(:,:,tt)/sum(bkReal(:,:,tt));
            end

            fkReal = ones(sqrt(M),1,J); %f0
            for kk = 1 : J

                gammak = 1/sqrt(M)*ones(sqrt(M));
                Pk=gammak.*exp(-((RealyBN(J-kk+1)-trellisOutput(:,:,J-kk+1))).^2/noiseVar);
                LckytmpBN(1,RealQm*(J-kk)+1) = log(fkReal(:,:,kk)'*(Afirst1.*Pk)*bkReal(:,:,kk+1)/(fkReal(:,:,kk)'*(Afirst0.*Pk)*bkReal(:,:,kk+1)));
                LckytmpBN(1,RealQm*(J-kk)+2) = log(fkReal(:,:,kk)'*(Asecond1.*Pk)*bkReal(:,:,kk+1)/(fkReal(:,:,kk)'*(Asecond0.*Pk)*bkReal(:,:,kk+1)));
                LckytmpBN(1,RealQm*(J-kk)+3) = log(fkReal(:,:,kk)'*(Athird1.*Pk)*bkReal(:,:,kk+1)/(fkReal(:,:,kk)'*(Athird0.*Pk)*bkReal(:,:,kk+1)));
                LckytmpBN(1,RealQm*(J-kk)+4) = log(fkReal(:,:,kk)'*(Aforth1.*Pk)*bkReal(:,:,kk+1)/(fkReal(:,:,kk)'*(Aforth0.*Pk)*bkReal(:,:,kk+1)));

                fkReal(:,:,kk+1) = Pk'*fkReal(:,:,kk);   % 前向更新 fk
                fkReal(:,:,kk+1) = fkReal(:,:,kk+1) / sum(fkReal(:,:,kk+1));
            end

            %% imag
            % Leimag = Leckp(2:2:end);
            bkImag =  ones(sqrt(M),1,J+1);  %b_Nt+1
            for tt = J:-1:2   %这里的Nt理解为状态图从后往前
                gammakpost = 1/sqrt(M)*ones(sqrt(M));
                PkPost = gammakpost.* (exp(-((ImagyBN(J-tt+1)-trellisOutput(:,:,J-tt+1))).^2/noiseVar)); % Pk,gamma构成的P矩阵

                bkImag(:,:,tt) = PkPost*bkImag(:,:,tt+1);
                bkImag(:,:,tt) = bkImag(:,:,tt)/sum(bkImag(:,:,tt));
            end
            fkImag = ones(sqrt(M),1,J); %f0
            for kk = 1 : J

                gammak = 1/sqrt(M)*ones(sqrt(M));

                Pk=gammak.*exp(-((ImagyBN(J-kk+1)-trellisOutput(:,:,J-kk+1))).^2/noiseVar);
                LckytmpBN(2,RealQm*(J-kk)+1) = log(fkImag(:,:,kk)'*(Afirst1.*Pk)*bkImag(:,:,kk+1)/(fkImag(:,:,kk)'*(Afirst0.*Pk)*bkImag(:,:,kk+1)));
                LckytmpBN(2,RealQm*(J-kk)+2) = log(fkImag(:,:,kk)'*(Asecond1.*Pk)*bkImag(:,:,kk+1)/(fkImag(:,:,kk)'*(Asecond0.*Pk)*bkImag(:,:,kk+1)));
                LckytmpBN(2,RealQm*(J-kk)+3) = log(fkImag(:,:,kk)'*(Athird1.*Pk)*bkImag(:,:,kk+1)/(fkImag(:,:,kk)'*(Athird0.*Pk)*bkImag(:,:,kk+1)));
                LckytmpBN(2,RealQm*(J-kk)+4) = log(fkImag(:,:,kk)'*(Aforth1.*Pk)*bkImag(:,:,kk+1)/(fkImag(:,:,kk)'*(Aforth0.*Pk)*bkImag(:,:,kk+1)));

                fkImag(:,:,kk+1) = Pk'*fkImag(:,:,kk);   % 前向更新 fk
                fkImag(:,:,kk+1) = fkImag(:,:,kk+1) / sum(fkImag(:,:,kk+1));
            end
            %%  test end
            Lckytmp(:,RealQm*J*(bn-1)+1:RealQm*J*bn) = LckytmpBN;
        end
        Lcky(:,nn) = reshape(Lckytmp,[],1);
        %% ================================== end

        % trellisOutput= ones(sqrt(M),sqrt(M),Nt);
        %% 计算状态转移矩阵  虚部实部共用一个，因为信道和状态都是给定的
        % for ii = Nt:-1:1
        %     if ii == Nt
        %         for mm = 1:sqrt(M)
        %             trellisOutput(:,mm,ii) = bidiagMat(ii,ii)*state(mm); %第Nt根天线的状态图输出
        %         end
        %     else
        %         for mm = 1:sqrt(M)
        %             trellisOutput(:,mm,ii) = bidiagMat(ii,ii+1)*state + bidiagMat(ii,ii)*state(mm); %第Nt根天线的状态图输出
        %         end
        %     end
        % end
        %% real
        % % 后向传递状态信息
        % bkReal =  ones(sqrt(M),1,Nt+1);  %b_Nt+1
        % LeReal = Leckp(1:2:end);
        % for tt = Nt:-1:2   %这里的Nt理解为状态图从后往前
        %     gammakpost = exp(Afirst1.*LeReal(Qm/2*(Nt-tt+1)-(Qm/2-1)))/(1 + exp(LeReal(Qm/2*(Nt-tt+1)-(Qm/2-1)))).*...
        %         (exp(Asecond1.*LeReal(Qm/2*(Nt-tt+1)-(Qm/2-2)))/(1 + exp(LeReal(Qm/2*(Nt-tt+1)-(Qm/2-2))))).*...
        %         (exp(Athird1.*LeReal(Qm/2*(Nt-tt+1)-(Qm/2-3)))/(1 + exp(LeReal(Qm/2*(Nt-tt+1)-(Qm/2-3))))).*...
        %         (exp(Aforth1.*LeReal(Qm/2*(Nt-tt+1)))/(1 + exp(LeReal(Qm/2*(Nt-tt+1)))));
        %     PkPost = gammakpost.* (exp(-((Realy(Nt-tt+1)-trellisOutput(:,:,Nt-tt+1))).^2/noiseVar)); % Pk,gamma构成的P矩阵
        % 
        %     bkReal(:,:,tt) = PkPost*bkReal(:,:,tt+1);
        %     bkReal(:,:,tt) = bkReal(:,:,tt)/sum(bkReal(:,:,tt));
        % end
        % fkReal = ones(sqrt(M),1,Nt); %f0
        % for kk = 1 : Nt
        % 
        %     gammak = exp(Afirst1.*LeReal(Qm/2*(Nt-kk+1)-(Qm/2-1)))/(1 + exp(LeReal(Qm/2*(Nt-kk+1)-(Qm/2-1)))).*...
        %         exp(Asecond1.*LeReal(Qm/2*(Nt-kk+1)-(Qm/2-2)))/(1 + exp(LeReal(Qm/2*(Nt-kk+1)-(Qm/2-2)))).*...
        %         exp(Athird1.*LeReal(Qm/2*(Nt-kk+1)-(Qm/2-3)))/(1 + exp(LeReal(Qm/2*(Nt-kk+1)-(Qm/2-3)))).*...
        %         exp(Aforth1.*LeReal(Qm/2*(Nt-kk+1)))/(1 + exp(LeReal(Qm/2*(Nt-kk+1))));
        % 
        %     Pk=gammak.*exp(-((Realy(Nt-kk+1)-trellisOutput(:,:,Nt-kk+1))).^2/noiseVar);
        %     Lckytmp(1,RealQm*(Nt-kk)+1) = log(fkReal(:,:,kk)'*(Afirst1.*Pk)*bkReal(:,:,kk+1)/(fkReal(:,:,kk)'*(Afirst0.*Pk)*bkReal(:,:,kk+1)));
        %     Lckytmp(1,RealQm*(Nt-kk)+2) = log(fkReal(:,:,kk)'*(Asecond1.*Pk)*bkReal(:,:,kk+1)/(fkReal(:,:,kk)'*(Asecond0.*Pk)*bkReal(:,:,kk+1)));
        %     Lckytmp(1,RealQm*(Nt-kk)+3) = log(fkReal(:,:,kk)'*(Athird1.*Pk)*bkReal(:,:,kk+1)/(fkReal(:,:,kk)'*(Athird0.*Pk)*bkReal(:,:,kk+1)));
        %     Lckytmp(1,RealQm*(Nt-kk)+4) = log(fkReal(:,:,kk)'*(Aforth1.*Pk)*bkReal(:,:,kk+1)/(fkReal(:,:,kk)'*(Aforth0.*Pk)*bkReal(:,:,kk+1)));
        % 
        % 
        %     fkReal(:,:,kk+1) = Pk'*fkReal(:,:,kk);   % 前向更新 fk
        %     fkReal(:,:,kk+1) = fkReal(:,:,kk+1) / sum(fkReal(:,:,kk+1));
        % end
        % 
        % %% imag
        % Leimag = Leckp(2:2:end);
        % bkImag =  ones(sqrt(M),1,Nt+1);  %b_Nt+1
        % for tt = Nt:-1:2   %这里的Nt理解为状态图从后往前
        %     gammakpost = exp(Afirst1.*Leimag(Qm/2*(Nt-tt+1)-(Qm/2-1)))/(1 + exp(Leimag(Qm/2*(Nt-tt+1)-(Qm/2-1)))).*...
        %         (exp(Asecond1.*Leimag(Qm/2*(Nt-tt+1)-(Qm/2-2)))/(1 + exp(Leimag(Qm/2*(Nt-tt+1)-(Qm/2-2))))).*...
        %         (exp(Athird1.*Leimag(Qm/2*(Nt-tt+1)-(Qm/2-3)))/(1 + exp(Leimag(Qm/2*(Nt-tt+1)-(Qm/2-3))))).*...
        %         (exp(Aforth1.*Leimag(Qm/2*(Nt-tt+1)))/(1 + exp(Leimag(Qm/2*(Nt-tt+1)))));
        %     PkPost = gammakpost.* (exp(-((Imagy(Nt-tt+1)-trellisOutput(:,:,Nt-tt+1))).^2/noiseVar)); % Pk,gamma构成的P矩阵
        % 
        %     bkImag(:,:,tt) = PkPost*bkImag(:,:,tt+1);
        %     bkImag(:,:,tt) = bkImag(:,:,tt)/sum(bkImag(:,:,tt));
        % end
        % fkImag = ones(sqrt(M),1,Nt); %f0
        % for kk = 1 : Nt
        % 
        %     gammak = exp(Afirst1.*Leimag(Qm/2*(Nt-kk+1)-(Qm/2-1)))/(1 + exp(Leimag(Qm/2*(Nt-kk+1)-(Qm/2-1)))).*...
        %         exp(Asecond1.*Leimag(Qm/2*(Nt-kk+1)-(Qm/2-2)))/(1 + exp(Leimag(Qm/2*(Nt-kk+1)-(Qm/2-2)))).*...
        %         exp(Athird1.*Leimag(Qm/2*(Nt-kk+1)-(Qm/2-3)))/(1 + exp(Leimag(Qm/2*(Nt-kk+1)-(Qm/2-3)))).*...
        %         exp(Aforth1.*Leimag(Qm/2*(Nt-kk+1)))/(1 + exp(Leimag(Qm/2*(Nt-kk+1))));
        % 
        %     Pk=gammak.*exp(-((Imagy(Nt-kk+1)-trellisOutput(:,:,Nt-kk+1))).^2/noiseVar);
        %     Lckytmp(2,RealQm*(Nt-kk)+1) = log(fkImag(:,:,kk)'*(Afirst1.*Pk)*bkImag(:,:,kk+1)/(fkImag(:,:,kk)'*(Afirst0.*Pk)*bkImag(:,:,kk+1)));
        %     Lckytmp(2,RealQm*(Nt-kk)+2) = log(fkImag(:,:,kk)'*(Asecond1.*Pk)*bkImag(:,:,kk+1)/(fkImag(:,:,kk)'*(Asecond0.*Pk)*bkImag(:,:,kk+1)));
        %     Lckytmp(2,RealQm*(Nt-kk)+3) = log(fkImag(:,:,kk)'*(Athird1.*Pk)*bkImag(:,:,kk+1)/(fkImag(:,:,kk)'*(Athird0.*Pk)*bkImag(:,:,kk+1)));
        %     Lckytmp(2,RealQm*(Nt-kk)+4) = log(fkImag(:,:,kk)'*(Aforth1.*Pk)*bkImag(:,:,kk+1)/(fkImag(:,:,kk)'*(Aforth0.*Pk)*bkImag(:,:,kk+1)));
        % 
        % 
        %     fkImag(:,:,kk+1) = Pk'*fkImag(:,:,kk);   % 前向更新 fk
        %     fkImag(:,:,kk+1) = fkImag(:,:,kk+1) / sum(fkImag(:,:,kk+1));
        % end
        % %%  test end
        % Lcky(:,nn) = reshape(Lckytmp,[],1);
    end
    LLRext = reshape(Lcky, size(Leck));
else
    error('Wrong demodulation mode.');
end




end



