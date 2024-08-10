clear;
% clc;
close all
% M = 2; % bi=1 or -1
M = 4; % M-QAM
tau = 0;
trainBits = 512; % training bits length
testBits = 2048; % testing bits length

channelLen = 1;
antenNum = 6; % number of antenna in SIMO
L = 3; % ISI channel paths L = 3
seed = 3;
rng(seed);
h_rayleigh = sqrt(1/2) *(randn(antenNum,L,channelLen) + 1i*randn(antenNum,L,channelLen));
snrdB = 2:2:18;
% snrdB = 14;
% snrdB = 5;
% snrdB = 20;
nMonte = 1;
h = 0.5; % BT=0.5
beta = 0.9;
lr = 0.15;
BER = zeros(12,length(snrdB));
viterbi_BER = zeros(12,length(snrdB));
viterbi_SER = zeros(12,length(snrdB));
bitsPerSym = log2(M);
%% ---------------------------Rayleigh Channel---------------------------- %%
for Hidx = 1 : channelLen
    fprintf(['\n Loop is %d at ', datestr(now,'HH:MM'), '\n'], Hidx);
    % % aktrain = [randi([0 1],trainBits-3,1);zeros(3,1)];
    % aktrain = randi([0 1],trainBits*bitsPerSym,1);
    aktrainSymbol = randi( [0,M-1],1,trainBits );
    aktrain = reshape(de2bi(aktrainSymbol,bitsPerSym,"left-msb").',1,[]).';
    % aktraintmp = [1 1 0 1 1 0 0 0 1 0 0 0 1 1 0 1];
    % aktrain = [];
    % for kkk = 1 : trainBits/16
    %     aktrain = [aktrain,[1 1 0 1 1 0 0 0 1 0 0 0 1 1 0 1]];
    % end
    % aktrain = [aktrain,1].';
    % label = phaseCompute_ISI(2*aktrain-1);
    % aktraintmp = [0 1 1 0 1 0];
    % bktrain = aktrain; % encoder (could be convolutional code,etc.)

    % coded
    bktrain = FEC_Coding(aktrain,2);
    % % label = phaseCompute_ISI(2*bktrain-1);
    label = QAM_label_compute_L3(bktrain,M);
        % tabulate(label2)
    % xktrain = gfsk_modulation(bktrain,h,tau);
    xktrain = qammod(aktrain, M, 'bin',InputType='bit', UnitAveragePower=true);
    % xktrain = qammod(bktrain, M, 'bin',InputType='bit', UnitAveragePower=true);
    
    % uncoded
    bktrain2 = aktrain;
    % % label2 = phaseCompute_ISI(2*bktrain2-1);
    label2 = QAM_label_compute_L3(bktrain2,M);
    % % xktrain2 = gfsk_modulation(bktrain2,h,tau);
    xktrain2 = qammod(bktrain2, M, 'bin',InputType='bit', UnitAveragePower=true);
    
    channelH_Vec = h_rayleigh(:,:,Hidx);
    channelH = channelH_Vec(1,:);
    aktestSymbol = randi( [0,M-1],1,testBits );
    aktest = reshape(de2bi(aktestSymbol,bitsPerSym,"left-msb").',1,[]).';
    bktest = aktest; % encoder
    xktest = qammod(bktest, M, 'bin',InputType='bit', UnitAveragePower=true);
    xktest2 = qammod(aktest, M, 'bin',InputType='bit', UnitAveragePower=true);
                vktrainSIMO_M4 = zeros(length(xktrain),size(channelH_Vec,1));
            vktrainSIMO2_M4 = zeros(length(xktrain2),size(channelH_Vec,1));
            
            
            % ISI channel output
            % SIMO M =4
            for ii = 1 : size(channelH_Vec,1)
                % encoded
                vktrainSIMO_M4(:,ii) = filter(channelH_Vec(ii,:),1,xktrain);
            end
            for ii = 1 : size(channelH_Vec,1)
                % uncoded
                vktrainSIMO2_M4(:,ii) = filter(channelH_Vec(ii,:),1,xktrain2);
            end
            % SIMO M = 6
            vktrainSIMO_M2 = vktrainSIMO_M4(:,1:4); 
            vktrainSIMO2_M2 = vktrainSIMO2_M4(:,1:4);
            %
            vktrain = filter(channelH,1,xktrain); % ISI channel training pilot output in SISO, encoded
            vktrain2 = filter(channelH,1,xktrain2); % ISI channel training pilot output in SISO, uncoded
            vktest = filter(channelH,1,xktest); % ISI channel testing sequence output in SISO
            vktest2 = filter(channelH,1,xktest2);
            % M = 4
            vktestSIMO_M4 = zeros(length(xktest),size(channelH_Vec,1));
            vktestSIMO2_M4 = zeros(length(xktest2),size(channelH_Vec,1));
            % ISI channel output in SIMO
            for ii = 1 : size(channelH_Vec,1)
                chantmp = channelH_Vec(ii,:);
                vktestSIMO_M4(:,ii) = filter(chantmp,1,xktest);
            end
            % no convolutional code
            for ii = 1 : size(channelH_Vec,1)
                chantmp = channelH_Vec(ii,:);
                vktestSIMO2_M4(:,ii) = filter(chantmp,1,xktest2);
            end
            % M = 2
            vktestSIMO_M2 = vktestSIMO_M4(:,1:4);
            vktestSIMO2_M2 = vktestSIMO2_M4(:,1:4);
    %% -----------------------------------SNR------------------------------------- %%
    for jj = 1 : length(snrdB)
        rng(jj)
        fprintf(['\n Inner iteration is %d at ', datestr(now,'HH:MM'), '\n'], jj);

        %% ---------------Monte Carlo Simulation-----------------%%
        for nn = 1 : nMonte
            rng(nn);


            % Generate GassianNoise in SISO and SIMO (Training),encoded
            GassianNoiseTrain_SIMO_M4 = (randn(size(vktrainSIMO_M4)) + 1j*randn(size(vktrainSIMO_M4)))/sqrt(2)*db2mag(-snrdB(jj));
            % GassianNoiseTrain_SIMO_M2 = GassianNoiseTrain_SIMO_M4(:,1:2);    
            GassianNoiseTrain = GassianNoiseTrain_SIMO_M4(:,1);
            % Generate GassianNoise in SISO and SIMO (Testing), encoded
            GassianNoiseTest_SIMO_M4 = (randn(size(vktestSIMO_M4)) + 1j*randn(size(vktestSIMO_M4)))/sqrt(2)*db2mag(-snrdB(jj));
            % GassianNoiseTest_SIMO_M2 = GassianNoiseTest_SIMO_M4(:,1:2);
            GassianNoiseTest = GassianNoiseTest_SIMO_M4(:,1);

            yktrainSIMO_M4 = vktrainSIMO_M4 + GassianNoiseTrain_SIMO_M4;
            % yktrainSIMO_M2 = vktrainSIMO_M2 + GassianNoiseTrain_SIMO_M2;
            yktrain = vktrain + GassianNoiseTrain;

            % Generate GassianNoise in SISO and SIMO (Training), uncoded
            GassianNoiseTrain_SIMO2_M4 = (randn(size(vktrainSIMO2_M4)) + 1j*randn(size(vktrainSIMO2_M4)))/sqrt(2)*db2mag(-snrdB(jj));
            % GassianNoiseTrain_SIMO2_M2 = GassianNoiseTrain_SIMO2_M4(:,1:2);
            GassianNoiseTrain2 = GassianNoiseTrain_SIMO2_M4(:,1);

            yktrainSIMO2_M4 = vktrainSIMO2_M4 + GassianNoiseTrain_SIMO2_M4;
            yktrainSIMO2_M2 = yktrainSIMO2_M4(:,1:4);
            % yktrainSIMO2_M2 = vktrainSIMO2_M2 + GassianNoiseTrain_SIMO2_M2;
            yktrain2 = vktrain2 + GassianNoiseTrain2;

            yktest = vktest + GassianNoiseTest;
            yktestSIMO_M4 = vktestSIMO_M4 + GassianNoiseTest_SIMO_M4;
            yktestSIMO_M2 = yktestSIMO_M4(:,1:4);
            % yktestSIMO_M2 = vktestSIMO_M2 + GassianNoiseTest_SIMO_M2;

            yktest2 = vktest2 + GassianNoiseTest(1:testBits);
            %% encoded
            % one-bit quantized training pilot zk_train in SISO
            zktrain_real = sign(real(yktrain)); zktrain_real = zktrain_real ./ sqrt(2);
            zktrain_imag = sign(imag(yktrain)); zktrain_imag = zktrain_imag ./ sqrt(2);
            zktrain = zktrain_real + 1i*zktrain_imag;

            % one-bit quantized training pilot zk_train in SIMO
            % M = 4
            zktrainSIMO_real_M4 = sign(real(yktrainSIMO_M4)); zktrainSIMO_real_M4 = zktrainSIMO_real_M4 ./ sqrt(2);
            zktrainSIMO_imag_M4 = sign(imag(yktrainSIMO_M4)); zktrainSIMO_imag_M4 = zktrainSIMO_imag_M4 ./ sqrt(2);
            zktrainSIMO_M4 = zktrainSIMO_real_M4 + 1i*zktrainSIMO_imag_M4;
            
            
            % % M = 2
            % zktrainSIMO_real_M2 = sign(real(yktrainSIMO_M2)); zktrainSIMO_real_M2 = zktrainSIMO_real_M2 ./ sqrt(2);
            % zktrainSIMO_imag_M2 = sign(imag(yktrainSIMO_M2)); zktrainSIMO_imag_M2 = zktrainSIMO_imag_M2 ./ sqrt(2);
            % zktrainSIMO_M2 = zktrainSIMO_real_M2 + 1i*zktrainSIMO_imag_M2;
            zktrainSIMO_M2 = zktrainSIMO_M4(:,1:4);
            %% uncoded
            % one-bit quantized training pilot zk_train in SISO
            zktrain_real2 = sign(real(yktrain2)); zktrain_real2 = zktrain_real2 ./ sqrt(2);
            zktrain_imag2 = sign(imag(yktrain2)); zktrain_imag2 = zktrain_imag2 ./ sqrt(2);
            zktrain2 = zktrain_real2 + 1i*zktrain_imag2;

            % one-bit quantized training pilot zk_train in SIMO
            zktrainSIMO_real2_M4 = sign(real(yktrainSIMO2_M4)); zktrainSIMO_real2_M4 = zktrainSIMO_real2_M4 ./ sqrt(2);
            zktrainSIMO_imag2_M4 = sign(imag(yktrainSIMO2_M4)); zktrainSIMO_imag2_M4 = zktrainSIMO_imag2_M4 ./ sqrt(2);
            zktrainSIMO2_M4 = zktrainSIMO_real2_M4 + 1i*zktrainSIMO_imag2_M4;
            zktrainSIMO2_M2 = zktrainSIMO2_M4(:,1:4);

            % % save Cross_Entropy.mat zktrainSIMO_real_M4 zktrainSIMO_imag_M4 label2 trainBits
            % % BLMMSE estimator in SISO, encoded
            % estedH_BLMMSE = Bussgang_LMMSE_ISI_L3(xktrain.', zktrain(3:end).', 1);
            % % BLMMSE estimator in SISO, uncoded
            % estedH_BLMMSE2 = Bussgang_LMMSE_ISI_L3(xktrain2.', zktrain2(3:end).', 1);
            % % % BLS estimator in SISO
            % % estedH_BLS = Bussgang_LeastSquare_ISI(xktrain.', zktrain(2:end).', 1);
            % 
            % % BLS estimator in SIMO
            % % estedH_SIMO_BLS = Bussgang_LeastSquare_ISI(xktrain.', zktrainSIMO(2:end,:).', length(channelH_Vec(:,1)));
            % 
            % BLMMSE estimator in SIMO, encoded
            estedH_SIMO_BLMMSE_M4 = Bussgang_LMMSE_ISI_L3(xktrain.', zktrainSIMO_M4(3:end,:).', length(channelH_Vec(:,1)));
            % BLMMSE estimator in SIMO, uncoded
            estedH_SIMO_BLMMSE2_M4 = Bussgang_LMMSE_ISI_L3(xktrain2.', zktrainSIMO2_M4(3:end,:).', length(channelH_Vec(:,1)));
            % % M = 2
            estedH_SIMO_BLMMSE_M2 = estedH_SIMO_BLMMSE_M4(1:4,:);
            estedH_SIMO_BLMMSE2_M2 = estedH_SIMO_BLMMSE2_M4(1:4,:);
            % estedH_SIMO_BLMMSE_M2 = Bussgang_LMMSE_ISI(xktrain.', zktrainSIMO_M4(2:end,1:2).', length(channelH_Vec(1:2,1)));
            % % BLMMSE estimator in SIMO, encoded
            % estedH_SIMO_BLMMSE2_M2 = Bussgang_LMMSE_ISI(xktrain2.', zktrainSIMO2_M4(2:end,1:2).', length(channelH_Vec(1:2,1)));
            % 
            % Testing bits
            zktest_real_M4 = sign(real(yktestSIMO_M4)); zktest_real_M4 = zktest_real_M4 ./ sqrt(2);
            zktest_imag_M4 = sign(imag(yktestSIMO_M4)); zktest_imag_M4 = zktest_imag_M4 ./ sqrt(2);
            zktestSIMO_M4 = zktest_real_M4 + 1i.*zktest_imag_M4;
            zktest = zktestSIMO_M4(:,1);
            zktestSIMO_M2 = zktestSIMO_M4(:,1:4);
            % uncoded
            yktestSIMO2_M4 = vktestSIMO2_M4 + GassianNoiseTest_SIMO_M4(1:testBits,:);
            zktest_real2_M4 = sign(real(yktestSIMO2_M4)); zktest_real2_M4 = zktest_real2_M4 ./ sqrt(2);
            zktest_imag2_M4 = sign(imag(yktestSIMO2_M4)); zktest_imag2_M4 = zktest_imag2_M4 ./ sqrt(2);
            zktestSIMO2_M4 = zktest_real2_M4 + 1i.*zktest_imag2_M4;
            zktest2 = zktestSIMO2_M4(:,1);
            % M = 2
            yktestSIMO2_M2 = yktestSIMO2_M4(:,1:4);
            zktestSIMO2_M2 = zktestSIMO2_M4(:,1:4);

            %--------------SIMO--M=2------------%
            % Perfect CSI
            [DecodeBit5,DecodeSym5] = viterbi_one_bit_SIMO_4QAM_detection_L3( zktestSIMO2_M2(3:end,:).', channelH_Vec(1:4,:), db2pow(-snrdB(jj)) );
            viterbi_BER(5,jj) = mean( DecodeBit5(5:end).' ~= aktest(7:end) ) + viterbi_BER(5,jj);
            viterbi_SER(5,jj) = mean( DecodeSym5(1:end) ~= aktestSymbol( 2:end) ) + viterbi_SER(5,jj);
            % BLMMSE
            [DecodeBit6,DecodeSym6] = viterbi_one_bit_SIMO_4QAM_detection_L3( zktestSIMO2_M2(3:end,:).', estedH_SIMO_BLMMSE2_M2, 1);
            viterbi_BER(6,jj) = mean( DecodeBit6(5:end).' ~= aktest(7:end) ) + viterbi_BER(6,jj);
            viterbi_SER(6,jj) = mean( DecodeSym6(1:end) ~= aktestSymbol( 2:end) ) + viterbi_SER(6,jj);
            % Quasi-NN
            [Re_hMtr2_M2,Im_hMtr2_M2,loss4_M2] = forward_update_SIMO_ISI_Channel_4QAM_L3(zktrainSIMO_real2_M4(3:end,1:4), zktrainSIMO_imag2_M4(3:end,1:4), label2, beta, lr, 32, trainBits);
            [DecodeBit7,DecodeSym7] = QNN_viterbi_one_bit_SIMO_detection_4QAM_L3(zktestSIMO2_M2(3:end,:),Re_hMtr2_M2,Im_hMtr2_M2);
            viterbi_BER(7,jj) = mean(DecodeBit7(5:end).' ~= aktest( 7:end) ) + viterbi_BER(7,jj);
            viterbi_SER(7,jj) = mean( DecodeSym7(1:end) ~= aktestSymbol( 2:end) ) + viterbi_SER(7,jj);
        end
    end
end

% BERAvg = zeros(size(BER));
BERAvg = BER ./ ( nMonte.*channelLen) ;
BER_Viterbi_Avg = viterbi_BER ./ ( nMonte.*channelLen) ;

%% plot encoded
% plt = {'r+-','gx-','bo-','k-',...
%         'r<--','gs--','b*--','k--',...
%         'r^-.','gd-.','bp-.','k-.'};
plt = {'rs-','g+-','bx-','k-',...
        'rs--','g+--','bx--','k--',...
        'rs-.','g+-.','bx-.','k-.'};

%% plot uncoded
figure(10)
for ii = 5 : 7
    semilogy(snrdB, BER_Viterbi_Avg(ii,:), plt{ii},'LineWidth',1,'MarkerSize',6);
    hold on
end
xlabel('SNR (dB)');
% xticks(-10:4:22);xlim([-10 26])
ylabel('BER');
yticks([1e-5 1e-4 1e-3 1e-2 1e-1 1e0])
grid on

legend("Genie-aided method",...
    "Conventional method with estimated CSI",...
    "Proposed Quasi-NN method")
                                                          
title('uncoded,viterbi')