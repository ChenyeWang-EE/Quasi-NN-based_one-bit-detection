function [ distance ] = distance_branch( branch,sig)

len = length(branch);
distance = zeros(1,len);
for kk = 1:len
    % distance(kk) = abs(branch(kk)).^2-2*real(sig'*branch(kk));
    distance(kk) = norm(branch(kk)-sig,2).^2;
    % distance(kk) = ( abs(branch(kk)-sig) ).^2/noiseVar;
    % distance3(kk) = ( branch(kk) - sig )'*( branch(kk) - sig );
end

end
