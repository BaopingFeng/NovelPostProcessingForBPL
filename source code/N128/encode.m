function x = encode(u);
N = size(u,1); 
n = log2(N);
if n==1
    x = [mod(u(1)+u(2),2); u(2)];
    return;
else
    x1 = encode(mod(u(1:N/2,1)+u(N/2+1:N,1),2));
    x2 = encode(u(N/2+1:N,1));
    x = [x1; x2];
end
return
    