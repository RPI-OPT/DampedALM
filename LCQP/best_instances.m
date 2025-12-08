count = zeros(3,10);
for i = 1:3
    numgrad = zeros(5,1);
    for j = 1:10
        numgrad(1) = max(out1{i,j}.gradnum);
        numgrad(2) = max(out2{i,j}.gradnum);
        numgrad(3) = max(out3{i,j}.gradnum);
        numgrad(4) = max(out4{i,j}.gradnum);
        numgrad(5) = max(out5{i,j}.gradnum);
        [a,b]      = min(numgrad);
        count(i,j) = b;
    end   
end