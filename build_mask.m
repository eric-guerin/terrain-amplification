function [ mask ] = build_mask ( size )
% build a polnomial mask matrix of size size*size and returns it
% size should be an even integer

offset = 1.0-1.0/size;

mask = zeros(size,size);
radius = (size-1)*0.5;
for i=1:size
    for j=1:size
        x = (i-1-radius)/radius;
        y = (j-1-radius)/radius;
        val = 1-offset*(x*x+y*y);
        if (val<0)
            val = 0;
        end
        mask(i,j) = val*val;
    end
end
end

