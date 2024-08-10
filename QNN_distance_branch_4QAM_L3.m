function distance = QNN_distance_branch_4QAM_L3( testdata,Re_h,Im_h)

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
yprob = zeros(length(testdata),M);
for ii = 1 : length(testdata)
    % [yprob(ii,:)] = QNN_ISI_ChannelH(testdata(ii), Re_h,Im_h, xkMtr);
    [yprob(ii,:)] = QNN_ISI_ChannelH_4QAM_L3(testdata(ii), Re_h,Im_h, xkMtr);
end
yprob = yprob + eps;
distance = -log(yprob);

end

