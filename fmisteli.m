function f = fmisteli(x)

if length(x)~=2
    error('Function is 2D')
end

f=(x(1)-5)^2+(x(2)-3)^2;
