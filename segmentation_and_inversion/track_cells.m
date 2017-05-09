function [small_lat result]=track_cells(initial,lat,inittime,endtime)
%this function tracks which cells is what at each time step of lattice.
%result(j,i) is the result of the cell number in the i'th timestep of the
%j'th cell in the initial timestep inittime.  small_lat is the new lattice
%in the lattmin format to be taken in by lattmin_convert to be pared down.

sizeinit=length(initial);
result(1:length(initial),1:endtime-inittime)=0;
result(:,1)=initial(:)+1; %re-index properly.
small_lat(1)=makenewlat(lat(inittime),result(:,1));
for i=inittime:(endtime-1)
    clear pos celldist;
    pos(1:sizeinit,1:2)=0;
    celldist(1:sizeinit,1:2)=0;
    for j=1:sizeinit
        %this is the centroid position of the j'th cell at timestep i.
        pos(j,1)=centerx(lat(i).cells{result(j,(i-inittime+1))},lat(i));
        pos(j,2)=centery(lat(i).cells{result(j,(i-inittime+1))},lat(i));
        mindist=100;mindex=1;mindelta=[100 100];
        clear temp;clear delta;
        temp(1:2)=[0 0];
        for k=2:length(lat(i+1).cells)
            temp(1)=centerx(lat(i+1).cells{k},lat(i+1));
            temp(2)=centery(lat(i+1).cells{k},lat(i+1));
            tempdist=(pos(j,1)-temp(1))^2+(pos(j,2)-temp(2))^2;
            delta(1)=temp(1)-pos(j,1);delta(2)=temp(2)-pos(j,2);
            if tempdist < mindist
                mindist=tempdist;mindex=k;
                mindelta(1)=delta(1);mindelta(2)=delta(2);
            end
        end
        celldist(j,1)=mindist;
        celldist(j,2)=mindex;
        ccellnum=celldist(j,2)-1;
        if mindelta(1) > 20 || mindelta(2) > 20
            fprintf('cell %d from timestep %d is not well tracked!\n',ccellnum,i);
        end
    end
    result(:,i-inittime+2)=celldist(:,2)';
    small_lat(i-inittime+2)=makenewlat(lat(i+1),result(:,i-inittime+2));
end
result(:,:)=result(:,:)-1;


%here is where we create the new lattice.
function newlat=makenewlat(lat,init_cells)
%lat is the original lattice indexing, init_cells are the current cells
%that we are tracking in the timestep of interest, and newlat is the
%re-indexed lattice in lattmin format.  the results of this should be
%plugged into lattmin_convert() to create a format usable by alberto's
%code, which will be processed again by my c code into my .info format.

%this creates newlat.cells{} and counts the initial number of bonds
%as well as creates a list of them, excluding the new cell0 that we are
%creating (for later part!)
countbond=0;
for i=1:length(init_cells)
    nbneighb=length(lat.cells{init_cells(i)});
    for j=1:nbneighb
        countbond=countbond+1; %counts the bonds
        %this line reindexes bonds and cells.
        newlat.cells{i+1}(j)=countbond;
        %this is a temporary holder for newlat.bonds. this will be changed
        %because of vertex reindexing as well as addition of more bonds.
        newlat.bonds(countbond,:)=lat.bonds(lat.cells{init_cells(i)}(j),:);
        %create the tempbondlist so we can deal with vertices
        tempbondlist(countbond)=lat.cells{init_cells(i)}(j);
    end
end
num_init_bond=countbond;

%create new vertex list and vertex structure
vertlist=[];
for i=1:num_init_bond
    vert=lat.bonds(tempbondlist(i),1);
    flag=0;
    for j=1:length(vertlist)
        if vert==vertlist(j)
            flag=1;break;
        end
    end
    if flag==0
        vertlist(length(vertlist)+1)=vert;
        newlat.verts(length(vertlist),:)=lat.verts(vert,:);
    end
end

%re-index the vertices in newlat.bonds, and make sure the proper bonds now
%border our new cell0.
cell0bonds=0;
for i=1:length(newlat.bonds)
    for j=1:length(vertlist)
        if newlat.bonds(i,1)==vertlist(j)
            newlat.bonds(i,1)=j;
            flag1=1;break; %break is necessary to prevent double relabel
        end
    end
    for j=1:length(vertlist)
        if newlat.bonds(i,2)==vertlist(j)
            newlat.bonds(i,2)=j;
            flag2=1;break; %break is necessary to prevent double relabel.
        end
    end
    flag=0;
    for k=1:length(init_cells)
        if newlat.bonds(i,4)+1==init_cells(k)
            flag=1;break;
        end
    end
    if flag==0
        newlat.bonds(i,4)=0;
        cell0bonds=cell0bonds+1;
        latestcell0bond=i;
        %the above is used as the initial starting point of the next
        %part where we actually create cell0 and the extra bonds.
    end
    %here we relabel the cell references on the rest of the bonds.
    for k=1:length(init_cells)
        if newlat.bonds(i,4)+1==init_cells(k)
            newlat.bonds(i,4)=k;
        end
        if newlat.bonds(i,3)+1==init_cells(k)
            newlat.bonds(i,3)=k;
        end
    end
end


%create cell0 list of bonds and the new bond structures it refers to
oldlength=length(newlat.bonds);
newlat.bonds(num_init_bond+1,1)=newlat.bonds(latestcell0bond,2);
newlat.bonds(num_init_bond+1,2)=newlat.bonds(latestcell0bond,1);
newlat.bonds(num_init_bond+1,3)=newlat.bonds(latestcell0bond,4);
newlat.bonds(num_init_bond+1,4)=newlat.bonds(latestcell0bond,3);
cell0(1)=num_init_bond+1;
for i=2:cell0bonds
    for j=1:oldlength
        if (newlat.bonds(j,4)==0 && newlat.bonds(j,2)==newlat.bonds(latestcell0bond,1))
            latestcell0bond=j;break;
        end
    end
    cell0(i)=num_init_bond+i;
    newlat.bonds(num_init_bond+i,1)=newlat.bonds(latestcell0bond,2);
    newlat.bonds(num_init_bond+i,2)=newlat.bonds(latestcell0bond,1);
    newlat.bonds(num_init_bond+i,3)=newlat.bonds(latestcell0bond,4);
    newlat.bonds(num_init_bond+i,4)=newlat.bonds(latestcell0bond,3);

end
newlat.cells{1}=cell0;

function x=centerx(cell,lat)

x=0;
for i=1:length(cell)
    x=x+lat.verts(lat.bonds(cell(i),1),1)/length(cell);
end


function y=centery(cell,lat)

y=0;
for i=1:length(cell)
    y=y+lat.verts(lat.bonds(cell(i),1),2)/length(cell);
end