Years = string(24:34);
Ny = length(Years);
Years = reshape(repmat(Years, [3, 1]), [1, 3*Ny]);
FilePerYear = string(1:3);
EndFileName = strcat(Years, "_", repmat(FilePerYear, [1, Ny]));
Files = strcat("Phobos_", EndFileName, ".txt");
N = length(Files);
Time = [];
Ephs = [];
for i = 1: N
    [M, T] = txt2matrix(Files(i));
    
    if i ~= 1
        T(1) = [];
        T = T+T0;
        M(1, :) = [];
    end

    Ephs = [Ephs; M];
    Time = [Time, T];
    T0 = Time(end);
    
end

% Time(1) = [];
AU = 149597870.7;
Ephs(:, 1) = Ephs(:, 1)*AU;
Ephs(:, 6) = Ephs(:, 6)*pi/180;
Ephs(:, 6) = unwrap(Ephs(:, 6));
Ephs(:, 6) = Ephs(:, 6)*180/pi;

save('PhobosEphs.mat', 'Ephs', 'Time');

