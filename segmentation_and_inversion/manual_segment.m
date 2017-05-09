function manual_segment()

close all;

global NAME
global PICTURENAME
global TYPE
global OUTPUTNAME
global fig
global vertcount
global cellcount
global clonecount
NAME='archimed_test';TYPE='tif'; %for screenshot
%load WingNRG.mat;num=22;NAME=sprintf('WingNRG%d',num);TYPE='gif'; %for direct import
PICTURENAME=sprintf('%s.%s',NAME,TYPE);
OUTPUTNAME=sprintf('%s_conf',NAME);
fig=1;
vertcount=0;
cellcount=0;
clonecount=0;

global vert
global cell


if nargin < 3
    fig=1;
end

img=imread(PICTURENAME,TYPE);%for screenshot
img=img(:,:,1); %weird thing with gimp tifs.
%wimage=wingimage;
%img=wimage(:,:,num); %for direct import
%load WingLattice2.mat;draw_cells(Lattice(num),fig,2);hold on;


%segmentation and watershed stuff
image1=img;
Seg = struct('WSTil',[]);
h = fspecial('gaussian', 30,3);
im = imfilter(image1,h);
im2 = edge(image1,'log',2/1000);
seD = strel('disk',10);
se0 = strel('line',2,0);
se90 = strel('line',2,90);
im2 = imerode(bwareaopen(imfill(imdilate(im2,seD),'holes'),10000),seD);
im3 = im;
im3(~im2) = Inf;
BW = imregionalmin(im3);
L = watershed(imimposemin(image1,BW));
image2 = image1;
image2(L == 0) = 255;
s = regionprops(L);
Areas = cat(1,s.Area);
Big = find(Areas > 1000);
OK = find(Areas< 1000);
for i = 1 : length(Big)
    L(L == Big(i)) = 1;
end
Seg.WSTil = imclearborder(L);
image1=img;
L = Seg.WSTil;
image1(L == 0) = 256;
L(bwlabel(imfill(imdilate(im2bw(256-image1),[se90,se0]),'holes'))==0) = 1;
image1=img;
image1(L == 0) = 256;
h_image = image(image1,'ButtonDownFcn',@image_callback);
%set(h_image,'AlphaData',0.5);


h_figure=figure(fig);
set(h_figure,'KeyPressFcn',@keypress_callback);
h_panel = imscrollpanel(h_figure,h_image);
api=iptgetapi(h_panel);
api.setMagnification(2.0);
h_inf=impixelinfo;set(h_inf,'Position',[0 20 50 20]);



function image_callback(obj,event)

global vertcount
global cellcount
global vert
global clonecount

norm=25.0; %resizes to be something instead of pixels.
pointer=get(gca,'CurrentPoint');
type=get(gcf,'SelectionType');

if strcmp(type,'normal')==1;
    vertcount=vertcount+1;
    NUM=sprintf('%d',vertcount);
    text(pointer(2,1),pointer(2,2),NUM,'BackgroundColor',[.7 .9 .7]);
    vert(vertcount).pos=[pointer(2,1), pointer(2,2), 0]/norm;
end
if strcmp(type,'extend')==1;
    cellcount=cellcount+1;
    NUM=sprintf('%d',cellcount);
    text(pointer(2,1),pointer(2,2),NUM,'BackgroundColor',[0 1 0]);
end
if strcmp(type,'open')==1;
    cellcount=cellcount+1;
    clonecount=clonecount+1;
    NUM=sprintf('%d',cellcount);
    text(pointer(2,1),pointer(2,2),NUM,'BackgroundColor',[1 0 0]);
end


function keypress_callback(src,evnt)

global PICTURENAME
global TYPE
global OUTPUTNAME
global fig
global vertcount
global cellcount
global vert
global cell

if evnt.Character == 'Q' %end and process
    fprintf('For the next part, please enter neighbors in a counter-clockwise fashion!\n')
    fprintf('Of course, cell0 has the special property where it looks clockwise to us.\n')
    for i=1:vertcount
        %do neighboring verts
        QUESTION=sprintf('Neighboring verts/cells to vertex %d (2x3 matrix): ',i);
        tmp=input(QUESTION);
        if length(tmp)~= 3
            fprintf('needs to have 3 neighbors!');
        end
        vert(i).nvert=(tmp(1,1:3)-[1 1 1]); %subtract to get indexing right
        vert(i).ncell=tmp(2,1:3);
    end
    for i=1:cellcount+1
        QUESTION=sprintf('Neighboring verts to cell %d (vector): ',i-1);
        tmp=input(QUESTION);
        nbneighb=length(tmp);
        cell(i).nbneighb=nbneighb;
        cell(i).nvert=(tmp-repmat(1,1,nbneighb)); %subtract to get indexing right
    end
    
    %corrections
    while 1==1
        QUESTION=sprintf('Are there any corrections you would like to make? (y/n): ');
        tmp=input(QUESTION,'s');
        if tmp=='y'
            QUESTION=sprintf('Vertex/Cell/Deltas/Quit? (v/c/d/q): ');
            tmpy=input(QUESTION,'s');
            switch tmpy
                case 'v'
                    QUESTION=sprintf('Which vertex #?: ');
                    tmpx=input(QUESTION);
                    QUESTION=sprintf('Neighboring verts/cells to vertex %d (2x3 matrix): ',tmpx);
                    temp=input(QUESTION);
                    if length(temp)~=3
                        fprintf('needs to have 3 neighbors!\n');
                    end
                    vert(tmpx).nvert=(temp(1,1:3)-[1 1 1]);
                    vert(tmpx).ncell=temp(2,1:3);
                case 'c'
                    QUESTION=sprintf('Which cell #?: ');
                    tmpx=input(QUESTION);
                    QUESTION=sprintf('Neighboring verts to cell %d (vector): ',tmpx);
                    temp=input(QUESTION);
                    nbneighb=length(temp);
                    cell(tmpx+1).nbneighb=nbneighb;
                    cell(tmpx+1).nvert=(temp-repmat(1,1,nbneighb));%get indexing right
                case 'd'
                    QUESTION=sprintf('Which cells are delta cells? ');
                    tmpx=input(QUESTION);
                    cell(tmpx(:)+1).clone=1;
                case 'q'
                    tmp='n';
            end
        end
        if tmp=='n'
            break;
        end
    end
    
    %writing stuff to alberto's conf file format here...
    fid=fopen(OUTPUTNAME,'w');
    fprintf(fid,'%d\n',vertcount);fprintf(fid,'%d\n',cellcount+1);
    for i=1:vertcount
        for j=1:3
            vert(i).dist(j)=sqrt((vert(i).pos(1)-vert(vert(i).nvert(j)+1).pos(1))^2 + ...
                (vert(i).pos(2)-vert(vert(i).nvert(j)+1).pos(2))^2);
        end
        fprintf(fid,'%f %f %d %f %f %f \n',vert(i).pos(1),vert(i).pos(2),vert(i).pos(3), ...
            vert(i).dist(1),vert(i).dist(2),vert(i).dist(3));
        fprintf(fid,'%d %d %d\n',vert(i).nvert(1),vert(i).nvert(2),vert(i).nvert(3));
        fprintf(fid,'%d %d %d\n',vert(i).ncell(1),vert(i).ncell(2),vert(i).ncell(3));
    end
    for i=1:cellcount+1
        nbneighb=cell(i).nbneighb;
        fprintf(fid,'%d\n',nbneighb);
        for j=1:nbneighb
            fprintf(fid,'%d\n',cell(i).nvert(j));
        end
    end
    fclose(fid);
            
end
if evnt.Character == 'R' %reset everything
    
    vert=[];
    vertcount=0;
    cell=[];
    cellcount=0;
    close all;
    
    %img = imread(PICTURENAME,TYPE); %for screenshot
    wimage=wingimage;img=wimage(:,:,num);
    h_image = image(img,'ButtonDownFcn',@image_callback);
    h_figure=figure(fig);
    set(h_figure,'KeyPressFcn',@keypress_callback);
    h_panel = imscrollpanel(h_figure,h_image);
    api=iptgetapi(h_panel);
    api.setMagnification(1);
end
if evnt.Character =='E' %edit something
    while 1==1
        QUESTION=sprintf('Are there any corrections you would like to make? (y/n): ');
        tmp=input(QUESTION,'s');
        if tmp=='y'
            QUESTION=sprintf('Vertex/Cell/Deltas/Quit? (v/c/d/q): ');
            tmpy=input(QUESTION,'s');
            switch tmpy
                case 'v'
                    QUESTION=sprintf('Which vertex #?: ');
                    tmpx=input(QUESTION);
                    QUESTION=sprintf('Neighboring verts/cells to vertex %d (2x3 matrix): ',tmpx);
                    temp=input(QUESTION);
                    if length(temp)~=3
                        fprintf('needs to have 3 neighbors!\n');
                    end
                    vert(tmpx).nvert=(temp(1,1:3)-[1 1 1]);
                    vert(tmpx).ncell=temp(2,1:3);
                case 'c'
                    QUESTION=sprintf('Which cell #?: ');
                    tmpx=input(QUESTION);
                    QUESTION=sprintf('Neighboring verts to cell %d (vector): ',tmpx);
                    temp=input(QUESTION);
                    nbneighb=length(temp);
                    cell(tmpx+1).nbneighb=nbneighb;
                    cell(tmpx+1).nvert=(temp-repmat(1,1,nbneighb));%get indexing right
                case 'd'
                    QUESTION=sprintf('Which cells are delta cells? ');
                    tmpx=input(QUESTION);
                    cell(tmpx(:)+1).clone=1;
                case 'q'
                    tmp='n';
            end
        end
        if tmp=='n'
            break;
        end
    end
end
if evnt.Character=='P' %process and print data
    %writing stuff to alberto's conf file format here...
    fid=fopen(OUTPUTNAME,'w');
    fprintf(fid,'%d\n',vertcount);fprintf(fid,'%d\n',cellcount+1);
    for i=1:vertcount
        for j=1:3
            vert(i).dist(j)=sqrt((vert(i).pos(1)-vert(vert(i).nvert(j)+1).pos(1))^2 + ...
                (vert(i).pos(2)-vert(vert(i).nvert(j)+1).pos(2))^2);
        end
        fprintf(fid,'%f %f %d %f %f %f \n',vert(i).pos(1),vert(i).pos(2),vert(i).pos(3), ...
            vert(i).dist(1),vert(i).dist(2),vert(i).dist(3));
        fprintf(fid,'%d %d %d\n',vert(i).nvert(1),vert(i).nvert(2),vert(i).nvert(3));
        fprintf(fid,'%d %d %d\n',vert(i).ncell(1),vert(i).ncell(2),vert(i).ncell(3));
    end
    numclone=0;
    for i=1:cellcount+1
        nbneighb=cell(i).nbneighb;
        fprintf(fid,'%d\n',nbneighb);
        if(cell(i).clone==1)
            numclone=numclone+1;
        end
        for j=1:nbneighb
            fprintf(fid,'%d\n',cell(i).nvert(j));
        end
    end
    fprintf(fid,'%d\n',numclone);
    for i=1:cellcount+1
        if cell(i).clone==1
            fprintf(fid,'%d\n',i);
        end
    end
    fclose(fid);
end
if evnt.Character=='I' %grab neighbor data from somewhere else
    INPUTNAME=input('Name of conf-type file to grab neighbor data from: ','s');
    fid=fopen(INPUTNAME,'r');
    n_vert=fscanf(fid,'%d');
    n_cell=fscanf(fid,'%d');
    if n_vert~=vertcount
        fprintf('The number of vertices from this file do not match the segmentation!\n');
    end
    if n_cell~=cellcount
        fprintf('The number of cells from this file do not match the segmentation!\n');
    end
    QUESTION=sprintf('Do you wish to import neighbor information from %s? (y/n): ',INPUTNAME);
    cont=input(QUESTION,'s');
    if cont=='n'
        return;
    end
    for i=1:n_vert
        for j=1:6
            fscanf(fid,'%f',1); %get rid of position data
        end
        for j=1:3
            dummy=fscanf(fid,'%d',1);
            vert(i).nvert(j)=dummy;
        end
        for j=1:3
            dummy=fscanf(fid,'%d',1);
            vert(i).ncell(j)=dummy;
        end
    end
    for i=1:n_cell
        nbneighb=fscanf(fid,'%d',1);
        for j=1:nbneighb
            cell(i).nvert(j)=fscanf(fid,'%d',1);
        end
    end
end
