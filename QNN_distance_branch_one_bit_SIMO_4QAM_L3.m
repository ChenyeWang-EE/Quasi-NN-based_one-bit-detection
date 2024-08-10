function distance = QNN_distance_branch_one_bit_SIMO_4QAM_L3( testdata,Re_hMtr,Im_hMtr)

len = size(testdata,1);
% xkMtr = [ exp(1j*pi*(-0.25-0.5)), exp(1j*pi*(0.25-0.5)), exp(1j*pi*(-0.25+0.5)), exp(1j*pi*(0.25+0.5)), exp(1j*pi*(0-0.25)), exp(1j*pi*(0+0.25)), exp(1j*pi*(1-0.25)), exp(1j*pi*(1+0.25)),...
%          exp(1j*pi*(0.5-0.25)), exp(1j*pi*(0.5+0.25)), exp(1j*pi*(-0.25-0.5)), exp(1j*pi*(0.25-0.5)), exp(1j*pi*(1-0.25)), exp(1j*pi*(1+0.25)), exp(1j*pi*(0-0.25)), exp(1j*pi*(0+0.25));
%          exp(1j*pi*(0-0.25)), exp(1j*pi*(0-0.25)), exp(1j*pi*(0+0.25)), exp(1j*pi*(0+0.25)), exp(1j*pi*(0.5-0.25)), exp(1j*pi*(0.5-0.25)), exp(1j*pi*(0.5+0.25)), exp(1j*pi*(0.5+0.25)),...
%          exp(1j*pi*(1-0.25)), exp(1j*pi*(1-0.25)), exp(1j*pi*(1+0.25)), exp(1j*pi*(1+0.25)), exp(1j*pi*(-0.5-0.25)), exp(1j*pi*(-0.5-0.25)), exp(1j*pi*(-0.5+0.25)), exp(1j*pi*(-0.5+0.25))];
map= [-1 / sqrt(2) + 1i / sqrt(2);
      -1 / sqrt(2) - 1i / sqrt(2);
       1 / sqrt(2) + 1i / sqrt(2);
       1 / sqrt(2) - 1i / sqrt(2)];
xkMtr = [ repmat( [map(1),map(2),map(3), map(4)], 1, 16 );
          repelem( [ map(1), map(2), map(3), map(4)], 16);
          repmat( repelem( [ map(1), map(2), map(3), map(4)], 4), 1, 4)];
M = 64;
yprob = zeros(len,M);
for ii = 1 : len
    % [yprob(ii,:)] = QNN_SIMO_ISI_Channel(testdata(ii,:), Re_hMtr,Im_hMtr,xkMtr);
    [yprob(ii,:)] = QNN_SIMO_ISI_Channel_L3(testdata(ii,:), Re_hMtr,Im_hMtr,xkMtr);
end
yprob = yprob + eps;
distance = -log(yprob);

end