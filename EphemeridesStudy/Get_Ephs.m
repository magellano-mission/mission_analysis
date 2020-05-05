

Files = ["Phobos24_26.txt", "Phobos26_28.txt", "Phobos28_30.txt", "Phobos30_32.txt", "Phobos32_34.txt"];
N = length(Files);

Ephs = [];
Time = 0;
for i = 1: N
    
    T0 = Time(end);
    [M, T] = txt2matrix(Files(i));
    T = T + T0;
    Ephs = [Ephs; M];
    Time = [Time, T];
    
end

Time(1) = [];
AU = 149597870.7;
Ephs(:, 1) = Ephs(:, 1)*AU;

save('PhobosEphs.mat', 'Ephs', 'Time');

