function [S, s_norm, s_dir ] = lineofsight(Y1,Y2)
%S = time variation of line 
S = Y2 - Y1;
s_norm = zeros(size(S,1),1);
s_dir = s_norm;
for jj=1:size(S,1)
s_norm(jj) = norm(S(jj,:));
end
for jj=1:size(S,1)
s_dir(jj) = acosd(dot((S(jj,:)./s_norm(jj)),(S(1,:)./s_norm(1))));
end

end

