profile on
profile clear
clear
clc
%% parameters settings
rng(31);
Nt = 8;Nr = 8;
load(strcat('CDL-A_',num2str(Nr),'×',num2str(Nt),'.mat'));
M = 16;
if M==4
    ModType = 'QPSK';
elseif M==16
    ModType = '16QAM';
end

Qm = log2(M);
NumBits = 128;
snrdB = 32:2:42;
MonteCalo = 1000;20000;20000;20000;
channelError =0; 1; %0 无误差；1 有误差
v_nCurves = [...          % Curves
    1 ...
    ];
Protsname = strvcat(  ...
    strcat("Iterative",'CBD'));
sequence = randperm(2*NumBits);
BER = zeros(length(v_nCurves),length(snrdB));
IterNum = 3;
BER1= zeros(IterNum,length(snrdB));
BER2= zeros(1,length(snrdB));
numblock = 4;
%% simulation begins
for mm = 1:MonteCalo
    if mod(mm,100) == 0
        fprintf([ '\n', 'MM = %d ', datestr(now), '\n'], mm);
    end
    for ii = 1:length(snrdB)
        alpha = db2pow(-snrdB(ii))*Nr;
        ak = [round(rand(NumBits-3,1));zeros(3,1)];
        [CodedBits, trel] = FEC_Coding(ak);
        % appDec = comm.APPDecoder('TrellisStructure',trel, ...
        %                          'Algorithm','True APP', ...
        %                          'CodedBitLLROutputPort',true);
        InterBits = CodedBits(sequence);
        symbol = nrSymbolModulate(InterBits,ModType);
        channelH =squeeze(channeH(:,:,mm));

        % [U,S,V] = bidiagonal(channelH); % bi-diagonal matrix, H = U*S*V';
        [U,S,V] = IP_CBD(channelH,numblock);
        % [U,S,V] = svd(channelH);

        Heff = channelH*V;
        xk = reshape(symbol, [Nt,length(symbol)/Nt]);
        noise = (randn(size(xk)) + 1j*randn(size(xk)))/sqrt(2)*sqrt(alpha);
        yk = channelH*V*xk + noise;

        % [Q_p,Rb_p] = matmodify(Heff);
        [Q_p,Rb_p] = qr(Heff);

        if (v_nCurves(1)==1)
            yHat1 = Q_p'*yk;

            LBn = zeros(length(symbol)*Qm, 1);
            for nn = 1:IterNum
                LeCkP = LBn;
                LCkY = MIMO_block_optimal_detection(yHat1,LeCkP, alpha,Rb_p, M, Nt,1);
                LeCkY = LCkY - LeCkP;
                LeBkY = deinterleaving(LeCkY, sequence);
                LeBkY(LeBkY>20)=20;
                LeBkY(LeBkY<-20)=-20;
                [LAkP, LBkY] = turbo_decode(LeBkY, 1/2);

                LAk = LBkY - LeBkY;  %Lext(bk|p)
                LBn = LAk(sequence);

                BER1(nn,ii) = mean([LAkP;zeros(3, 1)]~= ak)+BER1(nn,ii);
            end
        end


        
    end
end
BER= BER/MonteCalo;
BER1 = BER1/MonteCalo;
BER2 = BER2/MonteCalo;
%% plot
fig1 = figure;
v_stLegend = [];
set(fig1, 'WindowStyle', 'docked');
v_stPlotType = strvcat('-rs', '-bo', '-mx', '-rv', '--bo','--rv','-ko');

for aa=1:IterNum
    v_stLegend = strvcat(v_stLegend,  strcat('J=', num2str(numblock), ',',' ',' iter ',num2str(aa)));
    semilogy(snrdB, BER1(aa,:), v_stPlotType(aa,:),'LineWidth',1,'MarkerSize',6);
    hold on;
end
xlabel('SNR [dB]');
ylabel('BER');
title(strcat(num2str(Nt),'×',num2str(Nr),', Modulation order = ',num2str(ModType)))
grid on;
% if (v_nCurves(1)==1)
%     semilogy(snrdB, BER1,'-gv','LineWidth',1,'MarkerSize',6);
%     semilogy(snrdB, BER2,'-cv','LineWidth',1,'MarkerSize',6);
% end
legend(v_stLegend,'Location','SouthWest');