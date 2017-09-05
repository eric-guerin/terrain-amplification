function [ t ] = terrain_dilate( terrain, s )
%

t = zeros(size(terrain)+s);
radius = floor(s/2);
for i=1:size(t,1)
    for j=1:size(t,2)
        I = i-radius;
        J = j-radius;
        if I<1
            I = 1;
        elseif I>size(terrain,1)
            I = size(terrain,1);
        end
        if J<1
            J = 1;
        elseif J>size(terrain,2)
            J = size(terrain,2);
        end
        t(i,j) = terrain(I,J);
    end
end

end

