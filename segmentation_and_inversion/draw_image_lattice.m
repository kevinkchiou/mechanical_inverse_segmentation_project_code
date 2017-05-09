function draw_image_lattice(lat,img,fig)
% draws cells according to sebastian's data structure. numcell is a vector
% containing all the cells we want to track.

if nargin < 3
    fig=1;
end

mode=2;
figure(fig);
close(fig);
figure(fig);

himg=imshow(img);hold on;

if mode==1
    numcell=length(lat.cells);
    for i=1:numcell
        nbneighb=length(lat.cells{i});
        for j=1:nbneighb
            bnum=lat.cells{i}(j); %the global index of current bond
            if lat.bonds(bnum,3)~=0 && lat.bonds(bnum,4)~=0
                vert1x=lat.verts(lat.bonds(bnum,1),1);
                vert1y=lat.verts(lat.bonds(bnum,1),2);
                vert2x=lat.verts(lat.bonds(bnum,2),1);
                vert2y=lat.verts(lat.bonds(bnum,2),2);
                vx=[vert1x,vert2x]';vy=[vert1y,vert2y]';
                line(vx,vy,'Color',[1 0 0]);
            end
        end
    end
elseif mode==2
    numcell=length(lat.cells);
    for i=1:numcell
        nbneighb=length(lat.cells{i});
        avgxpos=0;avgypos=0;
        for j=1:nbneighb
            bnum=lat.cells{i}(j); %the global index of current bond
            vert1x=lat.verts(lat.bonds(bnum,1),1);
            vert1y=lat.verts(lat.bonds(bnum,1),2);
            vert2x=lat.verts(lat.bonds(bnum,2),1);
            vert2y=lat.verts(lat.bonds(bnum,2),2);
            vx=[vert1x,vert2x]';vy=[vert1y,vert2y]';
            line(vx,vy,'Color',[1 0 0]);
            avgxpos=avgxpos+vert1x/nbneighb;
            avgypos=avgypos+vert1y/nbneighb;
        end
        avgxpos=avgxpos-10;
        CELLNUM=sprintf('%d',i-1);
        if i>1
            text(avgxpos,avgypos,CELLNUM,'BackgroundColor',[0 1 0]);
        end
    end
end

if nargin>3
    ncell=length(trackcells);
    for j=1:ncell
        nbneighb=length(lat.cells{trackcells(j)+1});
        avgxpos=0;avgypos=0;
        for i=1:nbneighb
            bnum=lat.cells{trackcells(j)+1}(i);
            if lat.bonds(bnum,3)~=0 && lat.bonds(bnum,4)~=0
                vert1x=lat.verts(lat.bonds(bnum,1),1);
                vert1y=lat.verts(lat.bonds(bnum,1),2);
                vert2x=lat.verts(lat.bonds(bnum,2),1);
                vert2y=lat.verts(lat.bonds(bnum,2),2);
                vx=[vert1x,vert2x]';vy=[vert1y,vert2y]';
                line(vx,vy,'Color',[1 0 0]);
                avgxpos=avgxpos+vert1x/nbneighb;
                avgypos=avgypos+vert1y/nbneighb;
            end
        end
    avgxpos=avgxpos-10;
    CELLNUM=sprintf('%d',j);
    text(avgxpos,avgypos,CELLNUM);
    end
end

hfig=figure(fig);
h_panel = imscrollpanel(hfig,himg);
api=iptgetapi(h_panel);
api.setMagnification(2.0);
hinf=impixelinfo;set(hinf,'Position',[0 20 50 20]);