function [  ] = terrain_super_resolution( factor, exemplarhr, inputterrain, output, masksize, offset_analysis, offset_synthesis)
% perform terrain superesolution based on an exemplar
% factor : the amplification factor
% exemplarhr : the high resolution exemplar filename
% inputterrain : the input terrain filename
% output : the output terrain filename
% masksize : the diameter of the patch on the low resolution input terrain
% and exemplar
% offset_analysis : offset of the patch extracted from the low resolution
% exemplar
% offset_synthesis : offset that will be used for the synthesis
%
% Usage example:
% terrain_super_resolution(4,'grandcanyonhr.png','sketchlr.png',16,8,8);
%
% Be aware that giving too big exemplar or too small masksize and/or offset
% can give memory issues

tic

% read the high resolution exemplar

Iexemplarhr = imread(exemplarhr);
if size(size(Iexemplarhr),2) == 3
    Iexemplarhr = double(Iexemplarhr(:,:,1))+255*double(Iexemplarhr(:,:,2))+255*255*double(Iexemplarhr(:,:,3));
    Iexemplarhr = Iexemplarhr(:,:)/65535.0;
elseif isa(Iexemplarhr(1,1),'uint16')
    Iexemplarhr = double(Iexemplarhr)/255.0;
else
    Iexemplarhr = double(Iexemplarhr);
end
    
% the low resolution exemplar is a resizing of the high resolution exemplar
Iexemplarlr = imresize(Iexemplarhr,1.0/factor); 

% read the input terrain
Iinputterrain = imread(inputterrain);
if isa(Iinputterrain(1,1),'uint16')
    Iinputterrain = double(Iinputterrain)/255.0;
end
Iinputterrain = double(Iinputterrain);


% optimize dictionnary for the exemplar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mask = build_mask(masksize);
offset = offset_analysis;

d1 = floor((size(Iexemplarlr,1)-masksize)/offset);
d2 = floor((size(Iexemplarlr,2)-masksize)/offset);

% building mask matrix
Crec=zeros(size(Iexemplarlr));

for i=0:(d1-1)
    for j=0:(d2-1)
        Crec(i*offset+1:i*offset+masksize,j*offset+1:j*offset+masksize) = ...
            Crec(i*offset+1:i*offset+masksize,j*offset+1:j*offset+masksize) + mask;
    end
end


% building mean coefficients
exemplarmeans = zeros(d1,d2);

for i=0:(d1-1)
    for j=0:(d2-1)
        v=Iexemplarlr(i*offset+1:i*offset+masksize,j*offset+1:j*offset+masksize);
        exemplarmeans(i+1,j+1) = mean(v(:));
    end
end


% building terrain - mean
nl = 1;
X=zeros((d1)*(d2),masksize*masksize);

for i=0:(d1-1)
    for j=0:(d2-1)
        v=Iexemplarlr(i*offset+1:i*offset+masksize,j*offset+1:j*offset+masksize) ...
            - exemplarmeans(i+1,j+1)*ones(masksize,masksize);
        v = v.*mask;
        X(nl,:)=v(:)';
        nl = nl+1;
    end
end

X=double(X);

% the dictionnary is the whole number of atoms
Dinit = X';

% normalize atoms
nl = 1;
for i=1:size(Dinit,2)
    if norm(Dinit(:,i)) ~= 0
        D(:,nl) = Dinit(:,i)/norm(Dinit(:,i));
        nl = nl+1;
    end
end
n = nl-1;


% optimize terrain with the dictionnary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

offset = offset_synthesis;
Isyn = terrain_dilate(Iinputterrain,masksize*2);

d1 = floor((size(Isyn,1)-masksize)/offset);
d2 = floor((size(Isyn,2)-masksize)/offset);

% building mask matrix
Crec=zeros(size(Isyn));

for i=0:(d1-1)
    for j=0:(d2-1)
        Crec(i*offset+1:i*offset+masksize,j*offset+1:j*offset+masksize) = ...
            Crec(i*offset+1:i*offset+masksize,j*offset+1:j*offset+masksize) + mask;
    end
end

means = zeros(d1,d2);

Isketch = Isyn;

% building means
for i=0:(d1-1)
    for j=0:(d2-1)
        v=Isketch(i*offset+1:i*offset+masksize,j*offset+1:j*offset+masksize);
        means(i+1,j+1) = mean(v(:));
    end
end

% building terrain - mean
nl = 1;

X=zeros((d1)*(d2),masksize*masksize);

for i=0:(d1-1)
    for j=0:(d2-1)
        v=Isketch(i*offset+1:i*offset+masksize,j*offset+1:j*offset+masksize) ...
            - means(i+1,j+1)*ones(masksize,masksize);
        v = v.*mask;
        X(nl,:)=v(:)';
        nl = nl+1;
    end
end

X=double(X);

newcoeffs = omp(D,X',D'*D,1); % orthogonal matching pursuit with a sparsity of 1

fprintf('Number of patches : %d\n',size(newcoeffs,2));

% build the high resolution dictionnary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

masksizehr = masksize*factor;

maskhr = build_mask(masksizehr);
offsethra = offset_analysis*factor;
offsethrs = offset_synthesis*factor;



% building terrain - mean
d1 = floor((size(Iexemplarlr,1)-masksize)/offset_analysis);
d2 = floor((size(Iexemplarlr,2)-masksize)/offset_analysis);

nl = 1;
X=zeros((d1)*(d2),masksizehr*masksizehr);

noisemeanshr = zeros(d1,d2);

for i=0:(d1-1)
    for j=0:(d2-1)
        v=Iexemplarhr(i*offsethra+1:i*offsethra+masksizehr,j*offsethra+1:j*offsethra+masksizehr);
        noisemeanshr(i+1,j+1) = mean(v(:));
    end
end

for i=0:(d1-1)
    for j=0:(d2-1)
        v=Iexemplarhr(i*offsethra+1:i*offsethra+masksizehr,j*offsethra+1:j*offsethra+masksizehr) ...
            - noisemeanshr(i+1,j+1)*ones(masksizehr,masksizehr);
        v = v.*maskhr;
        X(nl,:)=v(:)';
        nl = nl+1;
    end
end

X=double(X);

% the dictionnary is the whole number of atoms
Dinithr = X';

D = [ ];

% normalize atoms
nl = 1;
for i=1:size(Dinithr,2)
    if  norm(Dinit(:,i)) ~= 0 
        D(:,nl) = Dinithr(:,i)/norm(Dinit(:,i));
        nl = nl+1;
    end
end
n = nl-1;

%%%%%%%%%%%%%%

d1 = floor((size(Isyn,1)-masksize)/offset_synthesis);
d2 = floor((size(Isyn,2)-masksize)/offset_synthesis);



Y = newcoeffs'*D';

Isynhr = zeros(size(Iinputterrain,1)*factor+2*masksizehr,size(Iinputterrain,2)*factor+2*masksizehr);

for i=0:(d1-1)
    for j=0:(d2-1)
        Isynhr(i*offsethrs+1:i*offsethrs+masksizehr,j*offsethrs+1:j*offsethrs+masksizehr) = ...
            Isynhr(i*offsethrs+1:i*offsethrs+masksizehr,j*offsethrs+1:j*offsethrs+masksizehr) ...
             + reshape(Y(i*(d2)+j+1,:),masksizehr,masksizehr) + means(i+1,j+1).*maskhr;
    end
end

% building mask matrix
Crechr=zeros(size(Isynhr));

for i=0:(d1-1)
    for j=0:(d2-1)
        Crechr(i*offsethrs+1:i*offsethrs+masksizehr,j*offsethrs+1:j*offsethrs+masksizehr) = ...
            Crechr(i*offsethrs+1:i*offsethrs+masksizehr,j*offsethrs+1:j*offsethrs+masksizehr) + maskhr;
    end
end

ind=find(Crechr ~=0);
Isynhr(ind) = Isynhr(ind) ./ Crechr(ind);

Isynhr = Isynhr(masksizehr+1:masksizehr+size(Iinputterrain,1)*factor,masksizehr+1:masksizehr+size(Iinputterrain,2)*factor);

Isynhr = Isynhr/255;

toc

Isynhr = uint16(floor(Isynhr*65535));
imwrite(Isynhr,output);

end


