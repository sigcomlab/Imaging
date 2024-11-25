function A_new = changeelements(A_old, n, a)

m = size(A_old,1);

selection = logical(triu(ones(m),n)-triu(ones(m),n+1));

A_new = A_old;

A_new(selection) = a;

end