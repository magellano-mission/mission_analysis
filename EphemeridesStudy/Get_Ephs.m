

Years = string(24:34);
Ny = length(Years);
Years = reshape(repmat(Years, [3, 1]), [1, 33]);
FilePerYear = string(1:3);
EndFileName = strcat(Years, "_", repmat(FilePerYear, [1, Ny]));
Files = strcat("Phobos_", EndFileName, ".txt");
N = length(Files);

Ephs = [];
Time = 0;
for i = 1: N
    
    T0 = Time(end);
    [M, T] = txt2matrix(Files(i));
    M(end, :) = [];
    T(end) = [];
    T = T + T0;
    Ephs = [Ephs; M];
    Time = [Time, T];
    
end

Time(1) = [];
AU = 149597870.7;
Ephs(:, 1) = Ephs(:, 1)*AU;

save('PhobosEphs.mat', 'Ephs', 'Time');

