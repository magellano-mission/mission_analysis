function N = coverageNumber(lat, lon, t, Y, theta)

[xx,yy,zz,~,~,~] = cone(rad2deg(theta),[0 0 0; Y(1) Y(2) Y(3)], norm(Y));

x = xx(:,end);
y = yy(:,end);
z = zz(:,end);

r = [x, y, z];
tt = ones(length(x), 1) * t;

[lat_S, lon_S] = groundTrack(r, tt);
lon_S = real(lon_S);

% if lat < max(lat_S) && lat > min(lat_S) && lon < max(lon_S) && lon > min(lon_S) % check this
%     N = 1;
%         A = 1;
% else
%     N = 0;
%             B = 1;
% end
% 
% if min(lon_S) < -150 && max(lon_S) > 150
%     if lon > 0
%         if lat < max(lat_S) && lat > min(lat_S) && lon > min(lon_S(lon_S>0)) % check this
%             N = 1;
%                     C = 1;
%         else
%             N = 0;
%                     D = 1;
%         end
%     else
%         if lat < max(lat_S) && lat > min(lat_S) && lon < max(lon_S(lon_S<0)) % check this
%             N = 1;
%                     E = 1;
%         else
%             N = 0;
%                     F = 1;
%         end
%     end
% end

if min(lon_S) < -150 && max(lon_S) > 150
    p = 0;

    if lon >= 0
        lat_S = lat_S(lon_S >= 0);
        lon_S = lon_S(lon_S >= 0);
        x = min(lon_S);
        if x > 25
            p = 1;
        end
    else
        lat_S = lat_S(lon_S < 0);
        lon_S = lon_S(lon_S < 0);
        x = max(lon_S);
        if x < - 50
            p = 1;
        end
    end

    if p == 0
        A = 0;
        B = sign(lon) * 180;
        C = sign(Y(3)) * 90;
        D = C;
        if lon >= 0
            [lon_S, I] = sort(lon_S);
            lat_S = lat_S(I);
            lon_S = [A; lon_S; B];
            lat_S = [C; lat_S; D];
        else
            [lon_S, I] = sort(lon_S);
            lat_S = lat_S(I);
            lon_S = [B; lon_S; A];
            lat_S = [C; lat_S; D];
        end
    else
        A = sign(lon) * 180;
        B = A;
        C = max(lat_S);
        D = min(lat_S);
        if lon >= 0
            [lon_S, I] = sort(lon_S);
            lat_S = lat_S(I);
            lon_S = [lon_S; A; B];
            lat_S = [lat_S; C; D];
        else
            [lon_S, I] = sort(lon_S);
            lat_S = lat_S(I);
            lon_S = [B; lon_S; A];
            lat_S = [C; lat_S; D];
        end
    end
end

 N = inpolygon(lon, lat, lon_S, lat_S);
 
end