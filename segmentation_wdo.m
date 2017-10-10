%-------------------------------------------------------------------------
% Please refer to the following journal article in your research papers:
% Z. Bayraktar, M. Komurcu, J. A. Bossard and D. H. Werner, "The Wind 
% Driven Optimization Technique and its Application in Electromagnetics," 
% IEEE Transactions on Antennas and Propagation, Volume 61, Issue 5, 
% pages 2745 - 2757, May 2013.
%-------------------------------------------------------------------------
%%
clc;
close all;
clear all;
warning('off', 'all')
format long g;
f1=1; %Zoom level

for input=[1:80]             %Image number 1 to 80
    disp(num2str(input));
    name=num2str(input);
    input_name=['im(' name ')'];
    I1=((imread(['D:\Image_folder_path\' input_name '.jpg']))); %Image name is in the formar im(22).jpg
    I=I1(:,60:(size(I1,2)-10));
    I=uint8(I);
    orig=I;
    dis=I;


            %% Image Dilation

            X=I>20;
            h = fspecial('gaussian',[15 15],0.01); 
            X=imfilter(X,h);
            Y=imfill(double(X),'holes');
            Z=bwareaopen(Y,400);

            SE=strel('disk',2);
            U=imdilate(Z,SE);

            %% Canny Edge Detection
            FinalEdge=edge(U,'canny');

            %%
            [inmi,in]=min(sum(I'));
            xi = in;

            yi = ceil(size(I,1)/2);


            for k=round(xi):1:m
                l=round(yi);
                    if (FinalEdge(k,l)==1)
                    r2=k;
                    c2=l;
                    break;
                    end
            end
            %% LineImage
            MaskImage=zeros(m,n);
            MaskImage(10:m-10,10:n-10)=1;
            LineImage=immultiply(FinalEdge,MaskImage);
            %% Lower ROI
            LowerROI=I(r2-30:r2+30,1:n); %30 pixels abobe and below the near wall are taken as the ROI
            


im=LowerROI;
[his, bin]=imhist(im);

% User defined WDO parameters:
npop = 20;              % population size.
nVar =3 ;               % Dimension of the problem.
maxit = 100;            % Maximum number of iterations.
param.RT = 3;			% RT coefficient.
param.g = 0.2;			% gravitational constant.
param.alp = 0.4;		% constants in the update eq.
param.c = 0.4;			% coriolis effect.
maxV = 70;              % maximum allowed speed.
varMin =  1;			% Lower dimension boundary.
varMax= numel(his);		% Upper dimension boundary.
%---------------------------------------------------------------

% Initialize WDO population, position and velocity:
p=randi([varMin,varMax],[1 nVar npop]);

for i=1:npop
    pos(i,:) = p(:,:,i); % Randomize velocity:
end

vel=zeros(npop, nVar);

% Evaluate initial population:
for K=1:npop,
   	pres(K,:) = otsu(his, pos(K,:));
end
%----------------------------------------------------------------

% Finding best air parcel in the initial population :
[globalpres,indx] = max(pres);
globalpos = pos(indx,:);
maxpres(1) = max(pres);			% maximum pressure


% Rank the air parcels:
[sorted_pres, rank_ind] = sort(pres, 'descend');
% Sort the air parcels:
pos = pos(rank_ind,:);
keepglob(1) = globalpres;



iter = 1;   % iteration counter
for ij = 2:maxit,
    	% Update the velocity:
    	for i=1:npop
		% choose random dimensions:
		a = randperm(nVar);        			
		% choose velocity based on random dimension:
    		velot(i,:) = vel(i,a);				
        	vel(i,:) = (1-param.alp)*vel(i,:)-(param.g*pos(i,:))+ ...
				    abs(1-1/i)*((globalpos-pos(i,:)).*param.RT)+ ...
				    (param.c*velot(i,:)/i);
    	end
    
        	% Check velocity:
        	vel = min(vel, maxV);
        	vel = max(vel, -maxV);
		% Update air parcel positions:
    		pos = pos + vel;
        	pos = round(min(pos, varMax));
        	pos = max(pos, varMin); 
		% Evaluate population: (Pressure)
		for K=1:npop,
				pres(K,:) = otsu(his, pos(K,:));
		end

    	%----------------------------------------------------
    	% Finding best particle in population
    	[maxpres,indx] = max(pres);
    	maxpos = pos(indx,:);           	% min location for this iteration
    	%----------------------------------------------------
    	% Rank the air parcels:
    	[sorted_pres, rank_ind] = sort(pres);
    	% Sort the air parcels position, velocity and pressure:
    	pos = pos(rank_ind,:);
    	vel = vel(rank_ind,:);
    	pres = sorted_pres;  
    
    	% Updating the global best:
    	better = maxpres > globalpres;
    	if better
        		globalpres = maxpres;             % initialize global minimum
        		globalpos = maxpos;
   	end
	% Keep a record of the progress:
    	keepglob(ij) = globalpres;    	
end

pressure = transpose(keepglob);
a=sort(globalpos);
levels=a-[ones(1,nVar)];
imq=imquantize(im,levels);
imq=imq-1;
imf=uint8(imq*(255/nVar));
[l, k, m, wid]=findlimits3(imf);
[row, col]=size(imf);

% Display the IMT layer on the original image
    for i=35:col-35
       dis(r2-30+k(i), i)=255;
       dis(r2-30+l(i), i)=255;
    end
figure;
imshow(dis);

end
