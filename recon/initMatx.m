function A = initMatx(I,J,N)
% initialize struct with matrix entries depending on x.

A = [];
Nstr = num2str(N);

for i = 1:I
    for j = 1:J
        istr = num2str(i);
        jstr = num2str(j);
        
        eval(['A.m'  istr jstr ' = zeros(' Nstr ',' Nstr ',' Nstr ');'])
    end
end


end