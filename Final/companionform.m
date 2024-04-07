function F = companionform(C)
    n = size(C{1}, 1);
    p = size(C,2);
    A = horzcat(C{:});
    I = eye(n*(p-1));
    F = [A;
        [I zeros(n*(p-1), n)]];

end
