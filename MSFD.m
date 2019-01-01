function [B, V, time] = MSFD(A, l)

[cols, rows, values] = find(A');
[n, d] = size(A);
AA = zeros(0,d);
B = [];
id = 0;
j=0;
tic;
for i = 1: n
    AA = [AA; A(i,:)];
    nnzA = nnz(AA);
    if nnzA >= l*d || size(AA, 1) == d || i == n
        BB = SpShrk(rows(id +1: id+nnzA)-j, cols(id+1:id+nnzA), values(id+1:id+nnzA), n, d, l);
        if nnz(B) == 0
            B = BB;
        else
            [B, V] = DnShrk([B; BB], l);
        end
        id = nnzA + id;
        AA = zeros(0,d);
        j = i;
    end
end
time = toc;

        
        
        
        
        