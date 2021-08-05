%% Dr. Denis Tsygankov (2021)
%% Integrative Systems Biology Lab, Wallace H. Coulter Department of Biomedical Engineering, Georgia Institute of Technology and Emory University SOM 

%% If you use any part of this script, please cite:
%
% DJ Marston et al. "Correcting artifacts in ratiometric biosensor imaging; an improved approach for dividing noisy signals"
% Frontiers in Cell and Developmental Biology (2021), doi: 10.3389/fcell.2021.685825

%% This interactive script generates a number of images to illustrate the NCF method (just follow the prompts). 

% You can use the two examples of data provided with this scrip or your own data, but make sure that 
% you have three images (numerator channel, denominator channel, and cell mask) of the same dimensions. 

% Your tif-files may have multiple frames or just one. If the channel files have multiple frames but the cell mask file has only one frame, 
% this one mask will be used for any frame of the signal channels that you specified.

% If you use the included data for the first time, start with the default parameters (these are the parameters used for Figure 8 in the paper).
% Notice that this script uses an optimization which is based on the cost function defined by Equation 5 in the paper.  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

%%
[fileFRET,pathFRET,indxFRET] = uigetfile('*.tif','Select the Numerator Channel');
[fileCFP,pathCFP,indxCFP] = uigetfile('*.tif','Select the Denominator Channel');
[fileMASK,pathMASK,indxMASK] = uigetfile('*.tif','Select the Cell Mask');

if indxFRET==0 || indxCFP==0 || indxMASK==0
    h = errordlg('Numerator, Denominator, and Mask must be selected','File Error','modal');
    return;
end

info = imfinfo([pathFRET fileFRET]);
numFRET = numel(info);
info = imfinfo([pathCFP fileCFP]);
numCFP = numel(info);
info = imfinfo([pathMASK fileMASK]);
numMASK = numel(info);

if numFRET==1 && numCFP==1 && numMASK==1
    FRET = double(imread([pathFRET fileFRET],1));   
    CFP = double(imread([pathCFP fileCFP],1));   
    MASK = double(imread([pathMASK fileMASK],1));
else
    fr = [];
    while isempty(fr)
        prompt = {'Enter the frame number:'};
        dlgtitle = 'Frame Selection';
        dims = [1 30];
        definput = {'1'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        if isempty(answer)
            fr = [];
        else
            fr = str2num(answer{1});
 
            if numMASK>1 && fr>0 && fr<=min([numFRET,numCFP,numMASK])
                FRET = double(imread([pathFRET fileFRET],fr));   
                CFP = double(imread([pathCFP fileCFP],fr));   
                MASK = double(imread([pathMASK fileMASK],fr));
            elseif numMASK==1 && fr>0 && fr<=min(numFRET,numCFP)
                FRET = double(imread([pathFRET fileFRET],fr));   
                CFP = double(imread([pathCFP fileCFP],fr));   
                MASK = double(imread([pathMASK fileMASK],1));
            else
                uiwait(errordlg('Frame number is inconsistent with the data files','Parameter Error','modal'));
                fr = [];
            end
        end
    end
end
    
figure('Position',get(0,'screensize'));
s1 = subplot(1,3,1);
imagesc(FRET);
axis equal; axis ij; %axis off;
title(['Numerator channel. Frame # ' num2str(fr)]);
set(s1,'XLim',[0.5 size(FRET,2)+0.5],'YLim',[0.5 size(FRET,1)+0.5])

s2 = subplot(1,3,2);
imagesc(CFP);
axis equal; axis ij; %axis off;
title(['Denominator channel. Frame # ' num2str(fr)]);
set(s2,'XLim',[0.5 size(CFP,2)+0.5],'YLim',[0.5 size(CFP,1)+0.5])

s3 = subplot(1,3,3);
imagesc(MASK);
axis equal; axis ij; %axis off;
if numMASK == 1
    title('Cell mask');
else
    title(['Cell mask. Frame # ' num2str(fr)]);
end
set(s3,'XLim',[0.5 size(MASK,2)+0.5],'YLim',[0.5 size(MASK,1)+0.5])

%%
prompt = {'Top:','Bottom:','Left:','Right:'};
dlgtitle = 'Background box';
dims = [1 30];
definput = {'230','295','20','96'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

box = zeros(4,1);
for i = 1:4
    box(i) = str2num(answer{i});
end

hold(s1,'on');
plot(s1,[box(3) box(3)],[box(1) box(2)],'g','LineWidth',3); 
plot(s1,[box(4) box(4)],[box(1) box(2)],'g','LineWidth',3); 
plot(s1,[box(3) box(4)],[box(1) box(1)],'g','LineWidth',3);        
plot(s1,[box(3) box(4)],[box(2) box(2)],'g','LineWidth',3); 

hold(s2,'on');
plot(s2,[box(3) box(3)],[box(1) box(2)],'g','LineWidth',3); 
plot(s2,[box(4) box(4)],[box(1) box(2)],'g','LineWidth',3); 
plot(s2,[box(3) box(4)],[box(1) box(1)],'g','LineWidth',3);        
plot(s2,[box(3) box(4)],[box(2) box(2)],'g','LineWidth',3); 

boxF = FRET(box(1):box(2),box(3):box(4)); 
bgF = mean(boxF(:));
sbgF = std(boxF(:));

boxC = CFP(box(1):box(2),box(3):box(4)); 
bgC = mean(boxC(:));
sbgC = std(boxC(:));
    
R_box = (FRET - bgF)./(CFP - bgC);

%%
prompt = {'From:','To:','Off edge distance:','Line-scan position:'};
dlgtitle = 'Specify NCF range and sub-mask';
dims = [1 30];
definput = {'460','800','10','244'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

v_min = str2num(answer{1});
v_max = str2num(answer{2});
off_edge = str2num(answer{3});
Y = str2num(answer{4});

v = v_min:v_max;
cost = zeros(size(v));
SE = strel('disk', off_edge);
subMASK = imerode(MASK,SE);

for i = 1:length(v)

    R_ncf = (FRET - v(i))./CFP;
    tmp = (R_ncf(subMASK==1) - R_box(subMASK==1)).^2;
    cost(i) = sqrt(sum(tmp)/length(tmp));
end

[cmin,ind] = min(cost);
optNCF = v(ind); 

R_NCF = (FRET - optNCF)./CFP;

R_BOX = R_box;
R_BOX(~MASK) = mean(R_box(MASK==1));

M = sum(MASK,2);
y1 = find(M>0,1,'first');
y2 = find(M>0,1,'last');
%Y = round(y1/2+y2/2);


figure('Position',get(0,'screensize'));
colormap(gray);
s1 = subplot(2,3,1);
hold on;
plot(v,cost,'b','LineWidth',2);
plot(optNCF,cmin,'bo','MarkerFaceColor','r');
title(['Optimization result: NCF = ', num2str(optNCF)]);
set(s1,'Box','on','XLim',[v_min v_max]);

s2 = subplot(2,3,2);
hold on;
imagesc(MASK+subMASK);
plot([0.5 size(MASK,2)+0.5],[Y Y],'r');
axis equal; axis ij; %axis off;
title(['Cell region ' num2str(off_edge) ' pixels away from the edge']);
set(s2,'XLim',[0.5 size(MASK,2)+0.5],'YLim',[0.5 size(MASK,1)+0.5]);

s3 = subplot(2,3,3);
im = R_BOX;
im(~MASK) = min(R_BOX(MASK==1));
imagesc(im);
axis equal; axis ij; %axis off;
title('BOX ratio (masked)');
set(s3,'XLim',[0.5 size(MASK,2)+0.5],'YLim',[0.5 size(MASK,1)+0.5]);

s4 = subplot(2,3,4);
hold on;
plot(R_BOX(Y,:),'r');
plot(R_NCF(Y,:),'b');
title('Ratios for box (red) and NCF (blue) methods along the line'); 
set(s4,'Box','on','XLim',[0.5 size(MASK,2)+0.5]);

s5 = subplot(2,3,5);
im = R_NCF;
%t = max(R_NCF(MASK==1));
%b = min(R_NCF(MASK==1));
%im(im < b - 0.1*(t-b)) = b - 0.1*(t-b);
%im(im > t + 0.1*(t-b)) = t + 0.1*(t-b);
imagesc(im);
axis equal; axis ij; %axis off;
title('NCF ratio (whole image)');
set(s5,'XLim',[0.5 size(MASK,2)+0.5],'YLim',[0.5 size(MASK,1)+0.5]);

s6 = subplot(2,3,6);
im = R_NCF;
im(~MASK) = min(R_NCF(MASK==1));
imagesc(im);
axis equal; axis ij; %axis off;
title('NCF ratio (masked)');
set(s6,'XLim',[0.5 size(MASK,2)+0.5],'YLim',[0.5 size(MASK,1)+0.5]);




%%