function [ distance ] = distance_branch_one_bit_L2( zk,channelH,noiseVar)



yprob = zeros(8,8);
trellisOutput = statechangeoutput(channelH);
sigma = sqrt(noiseVar);
%yprob
% for ii = 1 : len
switch zk
    case (1+1i) ./ sqrt(2)
        yprob = qfunc(-real(trellisOutput)./sigma).*qfunc(-imag(trellisOutput)./sigma);
    case (1-1i) ./ sqrt(2)
        yprob = qfunc(-real(trellisOutput)./sigma).*( 1-qfunc(-imag(trellisOutput)./sigma) );
    case (-1-1i) ./ sqrt(2)
        yprob = ( 1-qfunc(-real(trellisOutput)./sigma) ).*( 1-qfunc(-imag(trellisOutput)./sigma) );
    case (-1+1i) ./ sqrt(2)
        yprob = ( 1-qfunc(-real(trellisOutput)./sigma) ).*qfunc(-imag(trellisOutput)./sigma);
end
% end
yprob = yprob + eps;
S = [0 0 0 0 0 0 1 2;
    0 0 3 4 0 0 0 0;
    5 6 0 0 0 0 0 0;
    0 0 0 0 7 8 0 0;
    0 0 9 10 0 0 0 0;
    0 0 0 0 0 0 11 12;
    0 0 0 0 13 14 0 0;
    15 16 0 0 0 0 0 0];
R = (1:1:16);
[~,pos]=ismember(R,S);
branch = yprob(pos);
% distance = zeros(1,len);
% for kk = 1:len
%     % distance(kk) = abs(branch(kk)).^2-2*real(sig'*branch(kk));
%     distance(kk) = norm(branch(kk)-sig,2).^2;
%     % distance(kk) = ( abs(branch(kk)-sig) ).^2/noiseVar;
%     % distance3(kk) = ( branch(kk) - sig )'*( branch(kk) - sig );
%     distance(kk)
% end
distance = -log(branch);

end
