clear;
% clc;
close all
M = 2; % bi=1 or -1
tau = 0;
trainBits = 256; % training bits length
testBits = 5e2; % testing bits length

channelLen = 10;
antenNum = 4; % number of antenna in SIMO
L = 2; % ISI channel paths L = 2
seed = 10;
rng(seed);
h_rayleigh = sqrt(1/2) *(randn(antenNum,L,channelLen) + 1i*randn(antenNum,L,channelLen));
% snrdB = -10:2:30;
snrdB = -10:2:26;
% snrdB = 10;
nMonte = 1;
h = 0.5; % BT=0.5
beta = 0.9;
lr = 0.15;
BER = zeros(12,length(snrdB));
viterbi_BER = zeros(12,length(snrdB));
iternum = 1;
turErrBit_Inf = zeros(length(snrdB), iternum);
turErrBit_SIMO_inf_M4 = zeros(length(snrdB), iternum);
turErrBit_SIMO_inf_M2 = zeros(length(snrdB), iternum);
% SISO
turErrBit_Onebit = zeros(length(snrdB), iternum);
turErrBit_Onebit_est = zeros(length(snrdB), iternum);
turErrBit_Onebit_QNN = zeros(length(snrdB), iternum);
% SIMO M = 4
turErrBit_SIMO_M4 = zeros(length(snrdB), iternum);
turErrBit_SIMO_est_M4 = zeros(length(snrdB), iternum);
turErrBit_SIMO_QNN_M4 = zeros(length(snrdB), iternum);
% SIMO M = 2
turErrBit_SIMO_M2 = zeros(length(snrdB), iternum);
turErrBit_SIMO_est_M2 = zeros(length(snrdB), iternum);
turErrBit_SIMO_QNN_M2 = zeros(length(snrdB), iternum);
% nMonte1 = 50;
% nMonte2 = 50;
%% ---------------------------Rayleigh Channel---------------------------- %%
for Hidx = 1 : channelLen
    fprintf(['\n Loop is %d at ', datestr(now,'HH:MM'), '\n'], Hidx);
    aktrain = [randi([0 1],trainBits-3,1);zeros(3,1)];
    % aktrain = zeros(trainBits,1);
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
    bktrain = FEC_Coding(aktrain,M);
    label = phaseCompute_ISI(2*bktrain-1);
    %     tabulate(label)
    xktrain = gfsk_modulation(bktrain,h,tau);
    % uncoded
    bktrain2 = aktrain;
    label2 = phaseCompute_ISI(2*bktrain2-1);
    xktrain2 = gfsk_modulation(bktrain2,h,tau);

    channelH_Vec = h_rayleigh(:,:,Hidx);
    channelH = channelH_Vec(1,:);
    %% -----------------------------------SNR------------------------------------- %%
    for jj = 1 : length(snrdB)
        fprintf(['\n Inner iteration is %d at ', datestr(now,'HH:MM'), '\n'], jj);

        %% ---------------Monte Carlo Simulation-----------------%%
        for nn = 1 : nMonte
            rng(nn);
            %test
            aktest = [randi([0 1],testBits-3,1);zeros(3,1)];
            bktest = FEC_Coding(aktest,M);
            % bktest = aktest; % encoder
            xktest = gfsk_modulation(bktest,h,tau); % encoded
            xktest2 = gfsk_modulation(aktest,h,tau); % uncoded
            % M = 4
            vktrainSIMO_M4 = zeros(trainBits*M,length(channelH_Vec));
            vktrainSIMO2_M4 = zeros(trainBits,length(channelH_Vec));
            
            
            % ISI channel output
            % SIMO M =4
            for ii = 1 : length(channelH_Vec)
                % encoded
                vktrainSIMO_M4(:,ii) = filter(channelH_Vec(ii,:),1,xktrain);
            end
            for ii = 1 : length(channelH_Vec)
                % uncoded
                vktrainSIMO2_M4(:,ii) = filter(channelH_Vec(ii,:),1,xktrain2);
            end
            % SIMO M = 2
            vktrainSIMO_M2 = vktrainSIMO_M4(:,1:2); 
            vktrainSIMO2_M2 = vktrainSIMO2_M4(:,1:2);
            %
            vktrain = filter(channelH,1,xktrain); % ISI channel training pilot output in SISO, encoded
            vktrain2 = filter(channelH,1,xktrain2); % ISI channel training pilot output in SISO, uncoded
            vktest = filter(channelH,1,xktest); % ISI channel testing sequence output in SISO
            vktest2 = filter(channelH,1,xktest2);
            % M = 4
            vktestSIMO_M4 = zeros(testBits*M,length(channelH_Vec));
            vktestSIMO2_M4 = zeros(testBits,length(channelH_Vec));
            % ISI channel output in SIMO
            for ii = 1 : length(channelH_Vec)
                chantmp = channelH_Vec(ii,:);
                vktestSIMO_M4(:,ii) = filter(chantmp,1,xktest);
            end
            % no convolutional code
            for ii = 1 : length(channelH_Vec)
                chantmp = channelH_Vec(ii,:);
                vktestSIMO2_M4(:,ii) = filter(chantmp,1,xktest2);
            end
            % M = 2
            vktestSIMO_M2 = vktestSIMO_M4(:,1:2);
            vktestSIMO2_M2 = vktestSIMO2_M4(:,1:2);

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
            yktrainSIMO2_M2 = yktrainSIMO2_M4(:,1:2);
            % yktrainSIMO2_M2 = vktrainSIMO2_M2 + GassianNoiseTrain_SIMO2_M2;
            yktrain2 = vktrain2 + GassianNoiseTrain2;

            yktest = vktest + GassianNoiseTest;
            yktestSIMO_M4 = vktestSIMO_M4 + GassianNoiseTest_SIMO_M4;
            yktestSIMO_M2 = yktestSIMO_M4(:,1:2);
            % yktestSIMO_M2 = vktestSIMO_M2 + GassianNoiseTest_SIMO_M2;

            yktest2 = vktest2 + GassianNoiseTest(1:500);
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
            zktrainSIMO_M2 = zktrainSIMO_M4(:,1:2);
            %% uncoded
            % one-bit quantized training pilot zk_train in SISO
            zktrain_real2 = sign(real(yktrain2)); zktrain_real2 = zktrain_real2 ./ sqrt(2);
            zktrain_imag2 = sign(imag(yktrain2)); zktrain_imag2 = zktrain_imag2 ./ sqrt(2);
            zktrain2 = zktrain_real2 + 1i*zktrain_imag2;

            % one-bit quantized training pilot zk_train in SIMO
            zktrainSIMO_real2_M4 = sign(real(yktrainSIMO2_M4)); zktrainSIMO_real2_M4 = zktrainSIMO_real2_M4 ./ sqrt(2);
            zktrainSIMO_imag2_M4 = sign(imag(yktrainSIMO2_M4)); zktrainSIMO_imag2_M4 = zktrainSIMO_imag2_M4 ./ sqrt(2);
            zktrainSIMO2_M4 = zktrainSIMO_real2_M4 + 1i*zktrainSIMO_imag2_M4;
            zktrainSIMO2_M2 = zktrainSIMO2_M4(:,1:2);

            % save Cross_Entropy.mat zktrainSIMO_real_M4 zktrainSIMO_imag_M4 label2 trainBits
            % BLMMSE estimator in SISO, encoded
            estedH_BLMMSE = Bussgang_LMMSE_ISI(xktrain.', zktrain(2:end).', 1);
            % BLMMSE estimator in SISO, uncoded
            estedH_BLMMSE2 = Bussgang_LMMSE_ISI(xktrain2.', zktrain2(2:end).', 1);
            % % BLS estimator in SISO
            % estedH_BLS = Bussgang_LeastSquare_ISI(xktrain.', zktrain(2:end).', 1);

            % BLS estimator in SIMO
            % estedH_SIMO_BLS = Bussgang_LeastSquare_ISI(xktrain.', zktrainSIMO(2:end,:).', length(channelH_Vec(:,1)));

            % BLMMSE estimator in SIMO, encoded
            estedH_SIMO_BLMMSE_M4 = Bussgang_LMMSE_ISI(xktrain.', zktrainSIMO_M4(2:end,:).', length(channelH_Vec(:,1)));
            % BLMMSE estimator in SIMO, encoded
            estedH_SIMO_BLMMSE2_M4 = Bussgang_LMMSE_ISI(xktrain2.', zktrainSIMO2_M4(2:end,:).', length(channelH_Vec(:,1)));
            % % M = 2
            estedH_SIMO_BLMMSE_M2 = estedH_SIMO_BLMMSE_M4(1:2,:);
            estedH_SIMO_BLMMSE2_M2 = estedH_SIMO_BLMMSE2_M4(1:2,:);
            % estedH_SIMO_BLMMSE_M2 = Bussgang_LMMSE_ISI(xktrain.', zktrainSIMO_M4(2:end,1:2).', length(channelH_Vec(1:2,1)));
            % % BLMMSE estimator in SIMO, encoded
            % estedH_SIMO_BLMMSE2_M2 = Bussgang_LMMSE_ISI(xktrain2.', zktrainSIMO2_M4(2:end,1:2).', length(channelH_Vec(1:2,1)));
            
            % Testing bits
            zktest_real_M4 = sign(real(yktestSIMO_M4)); zktest_real_M4 = zktest_real_M4 ./ sqrt(2);
            zktest_imag_M4 = sign(imag(yktestSIMO_M4)); zktest_imag_M4 = zktest_imag_M4 ./ sqrt(2);
            zktestSIMO_M4 = zktest_real_M4 + 1i.*zktest_imag_M4;
            zktest = zktestSIMO_M4(:,1);
            zktestSIMO_M2 = zktestSIMO_M4(:,1:2);
            % uncoded
            yktestSIMO2_M4 = vktestSIMO2_M4 + GassianNoiseTest_SIMO_M4(1:500,:);
            zktest_real2_M4 = sign(real(yktestSIMO2_M4)); zktest_real2_M4 = zktest_real2_M4 ./ sqrt(2);
            zktest_imag2_M4 = sign(imag(yktestSIMO2_M4)); zktest_imag2_M4 = zktest_imag2_M4 ./ sqrt(2);
            zktestSIMO2_M4 = zktest_real2_M4 + 1i.*zktest_imag2_M4;
            zktest2 = zktestSIMO2_M4(:,1);
            % M = 2
            yktestSIMO2_M2 = yktestSIMO2_M4(:,1:2);
            zktestSIMO2_M2 = zktestSIMO2_M4(:,1:2);
            %% ----------------------estimated information bits using BCJR Algo.-------------------%%
            % % % turErrBit_Inf = Turbo_iter(yktest, channelH, turErrBit_Inf, aktest, db2pow(-snrdB(jj)),jj,0);
            % % % %% Turbo SISO
            % % % % perfect CSI
            % % % turErrBit_Onebit = Turbo_iter(zktest, channelH, turErrBit_Onebit, aktest, db2pow(-snrdB(jj)),jj,1);
            % % % % BLMMSE
            % % % turErrBit_Onebit_est = Turbo_iter(zktest, estedH_BLMMSE, turErrBit_Onebit_est, aktest, 1,jj,1);
            % % % % Turbo equalization and detection for Quasi-NN
            % % % [Re_hVec,Im_hVec,loss] = forward_update_ISI_Channel(zktrain_real(2:end), zktrain_imag(2:end), label, beta, lr, 8,trainBits);
            % % % turErrBit_Onebit_QNN = Turbo_iter_QNN_ISI(zktest, Re_hVec, Im_hVec, turErrBit_Onebit_QNN, aktest, jj);
            % % % 
            % % % % Prior Probabilities LLR for BCJR Algo.
            % % % % LeBkP = zeros(length(zktrain), 1);
            LeBkPtest = zeros(length(zktest), 1); % encoded
            LeBkPtest2 = zeros(length(zktest2), 1);% uncoded

            %-----------------------------SISO---------------------------------%
            %     % perfect CSI
            akHat1 = turbo_sym_detection_one_bit_ISI2(zktest2(2:end), LeBkPtest2, db2pow(-snrdB(jj)), channelH, 0);
            % BLMMSE estimator
            akHat2 = turbo_sym_detection_one_bit_ISI2(zktest2(2:end), LeBkPtest2, db2pow(-snrdB(jj)), estedH_BLMMSE2, 0);
            % BLS estimator
            % akHat3 = turbo_sym_detection_one_bit_ISI2(zktest(2:end), LeBkPtest, 1, estedH_BLS);

            % QNN
            [Re_hVec2,Im_hVec2,loss3] = forward_update_ISI_Channel(zktrain_real2(2:end), zktrain_imag2(2:end), label2, beta, lr, 8,trainBits);
            akHat3 = QNN_detection_ISI_Channel(LeBkPtest2, zktest2(2:end), Re_hVec2,Im_hVec2, 0);
            % Ideal ADC
            akHat4 = turbo_sym_detection_InfADC_ISI(yktest2(2:end), LeBkPtest, db2pow(-snrdB(jj)), channelH, 0);
            
            %--------------------------BER of SISO---------------------------%
            BER(1,jj) = mean( akHat1 ~= aktest(log2(M)+2:end-1) ) + BER(1,jj);
            BER(2,jj) = mean( akHat2 ~= aktest(log2(M)+2:end-1) ) + BER(2,jj);
            BER(3,jj) = mean( akHat3 ~= aktest(log2(M)+2:end-1) ) + BER(3,jj);
            BER(4,jj) = mean( akHat4 ~= aktest(log2(M)+2:end-1) ) + BER(4,jj);
            %% Viterbi
            % Viterbi perfect CSI
            Decode_onebit1 = viterbi_one_bit_detection(zktest2(2:end), channelH,db2pow(-snrdB(jj)));
            viterbi_BER(1,jj) = mean(Decode_onebit1(2:end).' ~= aktest(3:end)) + viterbi_BER(1,jj);
            % Viterbi Bussgang LMMSE
            Decode_onebit2 = viterbi_one_bit_detection(zktest2(2:end), estedH_BLMMSE2, 1);
            viterbi_BER(2,jj) = mean(Decode_onebit2(2:end).' ~= aktest(3:end)) + viterbi_BER(2,jj);   
            % Viterbi Quasi-NN
            Decode_onebit3 = QNN_viterbi_one_bit_detection(zktest2(2:end),Re_hVec2,Im_hVec2);
            viterbi_BER(3,jj) = mean(Decode_onebit3(2:end).' ~= aktest(3:end)) + viterbi_BER(3,jj);
            % Viterbi inf ADC
            DecodeBit_onebit4  = Viterbi_ISI_L2( yktest2(2:end),channelH);
            % DecodeBit = viterbi_gfskDemod_2h( yktest2,channelH);
            viterbi_BER(4,jj) = mean(DecodeBit_onebit4(2:end).' ~= aktest(3:end)) + viterbi_BER(4,jj);
            % % % % if snrdB(jj) <= 12
            % % %     %% Turbo SIMO M = 2
            % % %     turErrBit_SIMO_inf_M2 = Turbo_iter_SIMO(yktestSIMO_M2, channelH_Vec(1:2,:), turErrBit_SIMO_inf_M2, aktest,db2pow(-snrdB(jj)), jj, 0);
            % % %     %perfect CSI
            % % %     turErrBit_SIMO_M2 = Turbo_iter_SIMO(zktestSIMO_M2, channelH_Vec(1:2,:), turErrBit_SIMO_M2, aktest,db2pow(-snrdB(jj)), jj,1);
            % % %     % BLMMSE
            % % %     turErrBit_SIMO_est_M2 = Turbo_iter_SIMO(zktestSIMO_M2, estedH_SIMO_BLMMSE_M2, turErrBit_SIMO_est_M2, aktest,1, jj,1);
            % % %     [Re_hMtr_M2,Im_hMtr_M2,loss2_M2] = forward_update_SIMO_ISI_Channel(zktrainSIMO_real_M4(2:end,1:2), zktrainSIMO_imag_M4(2:end,1:2), label, beta, lr, 16, trainBits);
            % % %     turErrBit_SIMO_QNN_M2 = Turbo_iter_QNN_ISI_SIMO(zktestSIMO_M2, Re_hMtr_M2, Im_hMtr_M2, turErrBit_SIMO_QNN_M2, aktest, jj);
            % % % % end
            % % % % if snrdB(jj) >= 6
            % % % if snrdB(jj) <= 12
            % % %     %% Turbo SIMO M = 4
            % % %     turErrBit_SIMO_inf_M4 = Turbo_iter_SIMO(yktestSIMO_M4, channelH_Vec, turErrBit_SIMO_inf_M4, aktest,db2pow(-snrdB(jj)), jj, 0);
            % % %     %perfect CSI
            % % %     turErrBit_SIMO_M4 = Turbo_iter_SIMO(zktestSIMO_M4, channelH_Vec, turErrBit_SIMO_M4, aktest,db2pow(-snrdB(jj)), jj,1);
            % % %     % BLMMSE
            % % %     turErrBit_SIMO_est_M4 = Turbo_iter_SIMO(zktestSIMO_M4, estedH_SIMO_BLMMSE_M4, turErrBit_SIMO_est_M4, aktest,1, jj,1);
            % % %     [Re_hMtr,Im_hMtr,loss2] = forward_update_SIMO_ISI_Channel(zktrainSIMO_real_M4(2:end,:), zktrainSIMO_imag_M4(2:end,:), label, beta, lr, 16, trainBits);
            % % %     turErrBit_SIMO_QNN_M4 = Turbo_iter_QNN_ISI_SIMO(zktestSIMO_M4, Re_hMtr, Im_hMtr, turErrBit_SIMO_QNN_M4, aktest, jj);
            % % % end
            %% -----------------------------SIMO M=2---------------------------------%
            % SIMO perfect CSI
            akHat5 = turbo_sym_detection_one_bit_SIMO_ISI2(zktestSIMO2_M2(2:end,:).', LeBkPtest2, db2pow(-snrdB(jj)), channelH_Vec(1:2,:), 0);
            %SIMO BLMMSE estimated CSI
            akHat6 = turbo_sym_detection_one_bit_SIMO_ISI2(zktestSIMO2_M2(2:end,:).', LeBkPtest2, db2pow(-snrdB(jj)), estedH_SIMO_BLMMSE2_M2, 0);
            % %SIMO BLS estimated CSI
            % akHat9 = turbo_sym_detection_one_bit_SIMO_ISI2(zktestSIMO(2:end,:).', LeBkPtest, 1, estedH_SIMO_BLS);
            % SIMO QNN v2
            [Re_hMtr2_M2,Im_hMtr2_M2,loss4_M2] = forward_update_SIMO_ISI_Channel(zktrainSIMO_real2_M4(2:end,1:2), zktrainSIMO_imag2_M4(2:end,1:2), label2, beta, lr, 16, trainBits);
            akHat7 = QNN_detection_SIMO_ISI_Channel(LeBkPtest2, zktestSIMO2_M2(2:end,:), Re_hMtr2_M2,Im_hMtr2_M2, 0);
            % SIMO Ideal ADC
            akHat8 = turbo_sym_detection_InfADC_SIMO_ISI(yktestSIMO2_M2(2:end,:).', LeBkPtest, db2pow(-snrdB(jj)), channelH_Vec(1:2,:),0);
            
            %% -----------------------------SIMO M=4---------------------------------%
            % SIMO perfect CSI
            akHat9 = turbo_sym_detection_one_bit_SIMO_ISI2(zktestSIMO2_M4(2:end,:).', LeBkPtest2, db2pow(-snrdB(jj)), channelH_Vec, 0);
            %SIMO BLMMSE estimated CSI
            akHat10 = turbo_sym_detection_one_bit_SIMO_ISI2(zktestSIMO2_M4(2:end,:).', LeBkPtest2, db2pow(-snrdB(jj)), estedH_SIMO_BLMMSE2_M4, 0);
            % %SIMO BLS estimated CSI
            % akHat9 = turbo_sym_detection_one_bit_SIMO_ISI2(zktestSIMO(2:end,:).', LeBkPtest, 1, estedH_SIMO_BLS);
            % SIMO QNN v2
            [Re_hMtr2_M4,Im_hMtr2_M4,loss4] = forward_update_SIMO_ISI_Channel(zktrainSIMO_real2_M4(2:end,:), zktrainSIMO_imag2_M4(2:end,:), label2, beta, lr, 16, trainBits);
            akHat11 = QNN_detection_SIMO_ISI_Channel(LeBkPtest2, zktestSIMO2_M4(2:end,:), Re_hMtr2_M4,Im_hMtr2_M4, 0);
            % SIMO Ideal ADC
            akHat12 = turbo_sym_detection_InfADC_SIMO_ISI(yktestSIMO2_M4(2:end,:).', LeBkPtest, db2pow(-snrdB(jj)), channelH_Vec,0);

            % figure(jj)
            % semilogy(1:50,loss2,'LineWidth',1)
            % hold on
            % mean(loss2)
            %% Viterbi SIMO
            % -----------------------------SIMO M=2---------------------------------%
            % perfect CSI
            DecodeBit_onebit5 = viterbi_one_bit_SIMO_detection(zktestSIMO2_M2(2:end,:).', channelH_Vec(1:2,:),db2pow(-snrdB(jj)));
            viterbi_BER(5,jj) = mean( DecodeBit_onebit5(2:end).' ~= aktest(3:end) ) + viterbi_BER(5,jj);
            % Bussgang LMMSE
            DecodeBit_onebit6 = viterbi_one_bit_SIMO_detection(zktestSIMO2_M2(2:end,:).', estedH_SIMO_BLMMSE2_M2, 1);
            viterbi_BER(6,jj) = mean( DecodeBit_onebit6(2:end).' ~= aktest(3:end) ) + viterbi_BER(6,jj);
            % Quasi-NN
            Decode_onebit7 = QNN_viterbi_one_bit_SIMO_detection(zktestSIMO2_M2(2:end,:),Re_hMtr2_M2,Im_hMtr2_M2);
            viterbi_BER(7,jj) = mean(Decode_onebit7(2:end).' ~= aktest(3:end)) + viterbi_BER(7,jj);
            % inf ADC
            Decode_onebit8 = Viterbi_ISI_SIMO_L2( yktestSIMO2_M2(2:end,:).',channelH_Vec(1:2,:));
            viterbi_BER(8,jj) = mean(Decode_onebit8(2:end).' ~= aktest(3:end)) + viterbi_BER(8,jj);
            % -----------------------------SIMO M=4---------------------------------%
            % perfect CSI
            DecodeBit_onebit9 = viterbi_one_bit_SIMO_detection(zktestSIMO2_M4(2:end,:).', channelH_Vec,db2pow(-snrdB(jj)));
            viterbi_BER(9,jj) = mean( DecodeBit_onebit9(2:end).' ~= aktest(3:end) ) + viterbi_BER(9,jj);
            % Bussgang LMMSE
            DecodeBit_onebit10 = viterbi_one_bit_SIMO_detection(zktestSIMO2_M4(2:end,:).', estedH_SIMO_BLMMSE2_M4, 1);
            viterbi_BER(10,jj) = mean( DecodeBit_onebit10(2:end).' ~= aktest(3:end) ) + viterbi_BER(10,jj);
            % Quasi-NN
            Decode_onebit11 = QNN_viterbi_one_bit_SIMO_detection(zktestSIMO2_M4(2:end,:),Re_hMtr2_M4,Im_hMtr2_M4);
            viterbi_BER(11,jj) = mean(Decode_onebit11(2:end).' ~= aktest(3:end)) + viterbi_BER(11,jj);
            % inf ADC
            Decode_onebit12 = Viterbi_ISI_SIMO_L2( yktestSIMO2_M4(2:end,:).',channelH_Vec );
            viterbi_BER(12,jj) = mean(Decode_onebit12(2:end).' ~= aktest(3:end)) + viterbi_BER(12,jj);
            %--------------------------BER of SIMO M = 2---------------------------%
            BER(5,jj) = mean( akHat5 ~= aktest(log2(M)+2:end-1) ) + BER(5,jj);
            BER(6,jj) = mean( akHat6 ~= aktest(log2(M)+2:end-1) ) + BER(6,jj);
            BER(7,jj) = mean( akHat7 ~= aktest(log2(M)+2:end-1) ) + BER(7,jj);
            BER(8,jj) = mean( akHat8 ~= aktest(log2(M)+2:end-1) ) + BER(8,jj);
            %--------------------------BER of SIMO M = 4---------------------------%
            BER(9,jj) = mean( akHat9 ~= aktest(log2(M)+2:end-1) ) + BER(9,jj);
            BER(10,jj) = mean( akHat10 ~= aktest(log2(M)+2:end-1) ) + BER(10,jj);
            BER(11,jj) = mean( akHat11 ~= aktest(log2(M)+2:end-1) ) + BER(11,jj);
            BER(12,jj) = mean( akHat12 ~= aktest(log2(M)+2:end-1) ) + BER(12,jj);
            % end
        end
    end
end

% BERAvg = zeros(size(BER));
BERAvg = BER ./ ( nMonte.*channelLen) ;
BER_Viterbi_Avg = viterbi_BER ./ ( nMonte.*channelLen) ;
% BERAvg_InfADC = turErrBit_Inf ./ ( nMonte*testBits*channelLen ) ; % Infinite ADC
% BERAvg_Onebit = turErrBit_Onebit ./ ( nMonte*testBits*channelLen ) ;
% BERAvg_Onebit_ested = turErrBit_Onebit_est ./ ( nMonte*testBits*channelLen ) ;
% BERAvg_Onebit_QNN = turErrBit_Onebit_QNN ./ ( nMonte*testBits*channelLen ) ;
% % turbo BER M = 2
% BERAvg_InfADC_SIMO_M2 = turErrBit_SIMO_inf_M2 ./ ( nMonte*testBits*channelLen ) ; % Infinite ADC
% BERAvg_SIMO_M2 = turErrBit_SIMO_M2 ./ ( nMonte*testBits*channelLen ) ;
% BERAvg_SIMO_ested_M2 = turErrBit_SIMO_est_M2 ./ ( nMonte*testBits*channelLen ) ;
% BERAvg_SIMO_QNN_M2 = turErrBit_SIMO_QNN_M2 ./ ( nMonte*testBits*channelLen ) ;
% % turbo BER M = 4
% BERAvg_InfADC_SIMO_M4 = turErrBit_SIMO_inf_M4 ./ ( nMonte*testBits*channelLen ) ; % Infinite ADC
% BERAvg_SIMO_M4 = turErrBit_SIMO_M4 ./ ( nMonte*testBits*channelLen ) ;
% BERAvg_SIMO_ested_M4 = turErrBit_SIMO_est_M4 ./ ( nMonte*testBits*channelLen ) ;
% BERAvg_SIMO_QNN_M4 = turErrBit_SIMO_QNN_M4 ./ ( nMonte*testBits*channelLen ) ;

%% plot encoded
% plt = {'r+-','gx-','bo-','k-',...
%         'r<--','gs--','b*--','k--',...
%         'r^-.','gd-.','bp-.','k-.'};
plt = {'rs-','g+-','bx-','k-',...
        'rs--','g+--','bx--','k--',...
        'rs-.','g+-.','bx-.','k-.'};

% % %% plot SISO
% % figure(3)
% % % M = 1
% % semilogy(snrdB, BERAvg_Onebit(:,1), plt{1},'LineWidth',1,'MarkerSize',6);
% % hold on
% % semilogy(snrdB, BERAvg_Onebit_ested(:,1), plt{2},'LineWidth',1,'MarkerSize',6);
% % hold on
% % semilogy(snrdB, BERAvg_Onebit_QNN(:,1), plt{3},'LineWidth',1,'MarkerSize',6);
% % hold on
% % semilogy(snrdB, BERAvg_InfADC.', plt{4},'LineWidth',1,'MarkerSize',6);
% % hold on
% % 
% % % M = 2
% % semilogy(snrdB, BERAvg_SIMO_M2(:,1), plt{5},'LineWidth',1,'MarkerSize',6);
% % hold on
% % semilogy(snrdB, BERAvg_SIMO_ested_M2(:,1), plt{6},'LineWidth',1,'MarkerSize',6);
% % hold on
% % semilogy(snrdB, BERAvg_SIMO_QNN_M2(:,1), plt{7},'LineWidth',1,'MarkerSize',6);
% % hold on
% % semilogy(snrdB, BERAvg_InfADC_SIMO_M2.', plt{8},'LineWidth',1,'MarkerSize',6);
% % 
% % % M = 4
% % semilogy(snrdB, BERAvg_SIMO_M4(:,1), plt{9},'LineWidth',1,'MarkerSize',6);
% % hold on
% % semilogy(snrdB, BERAvg_SIMO_ested_M4(:,1), plt{10},'LineWidth',1,'MarkerSize',6);
% % hold on
% % semilogy(snrdB, BERAvg_SIMO_QNN_M4(:,1), plt{11},'LineWidth',1,'MarkerSize',6);
% % hold on
% % semilogy(snrdB, BERAvg_InfADC_SIMO_M4.', plt{12},'LineWidth',1,'MarkerSize',6);
% legend("Genie-aided method, M=1",...
%     "Conventional method with estimated CSI, M=1",...
%     "Proposed Quasi-NN method, M=1",'Genie-aided method, ideal ADC, M=1',...
%     "Genie-aided method, M=2",...
%     "Conventional method with estimated CSI, M=2",...
%     "Proposed Quasi-NN method, M=2",'Genie-aided method, ideal ADC, M=2',...
%     "Genie-aided method, M=4",...
%     "Conventional method with estimated CSI, M=4",...
%     "Proposed Quasi-NN method, M=4",'Genie-aided method, ideal ADC, M=4',...
%                                                       'Location','SouthWest')
% % legend("Genie-aided method",...
% %     "Conventional method with estimated CSI",...
% %     "Proposed Quasi-NN method",'Genie-aided method, ideal ADC')
                                                          

% % xlabel('SNR (dB)');
% % xticks(-10:4:26);xlim([-10 26])
% % ylabel('BER');
% % yticks([1e-5 1e-4 1e-3 1e-2 1e-1 1e0])
% % grid on
% % title('Convolutional code')
%% plot uncoded
figure(4)
for ii = 1 : size(BER,1)
    semilogy(snrdB, BERAvg(ii,:), plt{ii},'LineWidth',1,'MarkerSize',6);
    hold on
end
xlabel('SNR (dB)');
xticks(-10:4:26);xlim([-10 26])
ylabel('BER');
yticks([1e-5 1e-4 1e-3 1e-2 1e-1 1e0])
grid on
% legend("Genie-aided method, M=1",...
%     "Conventional method with estimated CSI, M=1",...
%     "Proposed Quasi-NN method, M=1",'Genie-aided method, ideal ADC, M=1',...
%     "Genie-aided method, M=2",...
%     "Conventional method with estimated CSI, M=2",...
%     "Proposed Quasi-NN method, M=2",'Genie-aided method, ideal ADC, M=2',...
%     "Genie-aided method, M=4",...
%     "Conventional method with estimated CSI, M=4",...
%     "Proposed Quasi-NN method, M=4",'Genie-aided method, ideal ADC, M=4',...
%                                                       'Location','SouthWest')
legend("Genie-aided method",...
    "Conventional method with estimated CSI",...
    "Proposed Quasi-NN method",'Genie-aided method, ideal ADC')
                                                          
title('uncoded')

%% plot uncoded
figure(5)
for ii = 1 : size(BER,1)
    semilogy(snrdB, BER_Viterbi_Avg(ii,:), plt{ii},'LineWidth',1,'MarkerSize',6);
    hold on
end
xlabel('SNR (dB)');
xticks(-10:4:26);xlim([-10 26])
ylabel('BER');
yticks([1e-5 1e-4 1e-3 1e-2 1e-1 1e0])
grid on
% legend("Genie-aided method, M=1",...
%     "Conventional method with estimated CSI, M=1",...
%     "Proposed Quasi-NN method, M=1",'Genie-aided method, ideal ADC, M=1',...
%     "Genie-aided method, M=2",...
%     "Conventional method with estimated CSI, M=2",...
%     "Proposed Quasi-NN method, M=2",'Genie-aided method, ideal ADC, M=2',...
%     "Genie-aided method, M=4",...
%     "Conventional method with estimated CSI, M=4",...
%     "Proposed Quasi-NN method, M=4",'Genie-aided method, ideal ADC, M=4',...
%                                                       'Location','SouthWest')
legend("Genie-aided method",...
    "Conventional method with estimated CSI",...
    "Proposed Quasi-NN method",'Genie-aided method, ideal ADC')
                                                          
title('uncoded,viterbi')