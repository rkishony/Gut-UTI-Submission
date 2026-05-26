function TF = seteq(A, B)
% Check if two sets are equal

TF = isempty(setxor(A,B));
end