function distance = QNN_distance_branch_one_bit_SIMO( testdata,Re_hMtr,Im_hMtr)

len = size(testdata,1);
xkMtr = [ exp(1j*pi*(-0.25-0.5)), exp(1j*pi*(0.25-0.5)), exp(1j*pi*(-0.25+0.5)), exp(1j*pi*(0.25+0.5)), exp(1j*pi*(0-0.25)), exp(1j*pi*(0+0.25)), exp(1j*pi*(1-0.25)), exp(1j*pi*(1+0.25)),...
         exp(1j*pi*(0.5-0.25)), exp(1j*pi*(0.5+0.25)), exp(1j*pi*(-0.25-0.5)), exp(1j*pi*(0.25-0.5)), exp(1j*pi*(1-0.25)), exp(1j*pi*(1+0.25)), exp(1j*pi*(0-0.25)), exp(1j*pi*(0+0.25));
         exp(1j*pi*(0-0.25)), exp(1j*pi*(0-0.25)), exp(1j*pi*(0+0.25)), exp(1j*pi*(0+0.25)), exp(1j*pi*(0.5-0.25)), exp(1j*pi*(0.5-0.25)), exp(1j*pi*(0.5+0.25)), exp(1j*pi*(0.5+0.25)),...
         exp(1j*pi*(1-0.25)), exp(1j*pi*(1-0.25)), exp(1j*pi*(1+0.25)), exp(1j*pi*(1+0.25)), exp(1j*pi*(-0.5-0.25)), exp(1j*pi*(-0.5-0.25)), exp(1j*pi*(-0.5+0.25)), exp(1j*pi*(-0.5+0.25))];
M = 16;
yprob = zeros(len,M);
for ii = 1 : len
    % [yprob(ii,:)] = QNN_SIMO_ISI_Channel(testdata(ii,:), Re_hMtr,Im_hMtr,xkMtr);
    [yprob(ii,:)] = QNN_SIMO_ISI_Channel_2(testdata(ii,:), Re_hMtr,Im_hMtr,xkMtr);
end
yprob = yprob + eps;
distance = -log(yprob);

end