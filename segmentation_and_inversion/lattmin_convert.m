function [newlat oldlat vdat cdat]=lattmin_convert(lattice,t)
% This function converts data taken in by lattmin into a format that can be
% read by my own code and then build the *info files via the c routines
% already created.  Cheap man's shortcut for someone short on time (since I
% procrastinate.  So stop doing that!!).

FNAME=sprintf('~/sebastian/config_%0.2d',t);
fid=fopen(FNAME,'w');
[newlat oldlat]= pare_data(lattice,t); %makes vertices trivalent
lat=newlat;
numvert=length(lat.verts);
numcell=length(lat.cells);
fprintf(fid,'%d\n%d\n',numvert,numcell);
[vdat cdat]=create_data(lat); %rewrites data into vdat and cdat format

%the following writes the data into alberto's format so that i can use my
%c code to extract bond structure too instead of rewriting it all.
for i=1:length(vdat)
    if vdat(i).nvert(1)<0 || vdat(i).nvert(2)<0 || vdat(i).nvert(3)<0
        fprintf('nvert1 = %d, nvert2 = %d, nvert3 = %d, i = %d\n', ...
            vdat(i).nvert(1),vdat(i).nvert(2),vdat(i).nvert(3),i);
        error('bad!');
    end
    dist1=sqrt((vdat(i).xcoord-vdat(vdat(i).nvert(1)+1).xcoord)^2 + ...
        (vdat(i).ycoord - vdat(vdat(i).nvert(1)+1).ycoord)^2);
    dist2=sqrt((vdat(i).xcoord-vdat(vdat(i).nvert(2)+1).xcoord)^2 + ...
        (vdat(i).ycoord - vdat(vdat(i).nvert(2)+1).ycoord)^2);
    dist3=sqrt((vdat(i).xcoord-vdat(vdat(i).nvert(3)+1).xcoord)^2 + ...
        (vdat(i).ycoord - vdat(vdat(i).nvert(3)+1).ycoord)^2);
    fprintf(fid,'%d %d 0 %f %f %f\n', vdat(i).xcoord,vdat(i).ycoord,dist1,dist2,dist3);
    fprintf(fid,'%d %d %d\n',vdat(i).nvert(1),vdat(i).nvert(2),vdat(i).nvert(3));
    fprintf(fid,'%d %d %d\n',vdat(i).ncell(1),vdat(i).ncell(2),vdat(i).ncell(3));
end
for i=1:length(cdat)
    nbneighb=cdat(i).nbneighb;
    fprintf(fid,'%d\n',nbneighb);
    for j=1:nbneighb
        fprintf(fid,'%d\n',cdat(i).nvert(j));
    end
end
fclose(fid);


function [vdat cdat]= create_data(lat)

numvert=length(lat.verts);
numcell=length(lat.cells);
vdat=repmat(struct('index',-1,'ncell',[-1 -1 -1], ...
    'nvert',[-1 -1 -1],'xcoord',0,'ycoord',0,'zcoord',0),1,numvert);
cdat=repmat(struct('nbneighb',-1,'nvert',[-1 -1 -1]),1,numcell-1);

for i=1:numcell
    %don't need to include cell0, since it's taken care of in other
    %algorithms.  also cells only need to know their nvert to be successful
    %in our other algorithms, so the remainder of the data is left empty
    nbneighb=length(lat.cells{i});
    cdat(i).nbneighb=nbneighb;
    for j=1:nbneighb
        %this has funny order in order to give c-clockwise orientation
        cdat(i).nvert(j)=lat.bonds(lat.cells{i}(nbneighb-j+1),1);
        bond=lat.cells{i}(j);
        fbvert=lat.bonds(bond,1);
        nbvert=lat.bonds(bond,2);
        %now determine the neighboring cell
        if lat.bonds(bond,4)==i-1
            error('something is wrong with the structure!');
        end
        ocell=lat.bonds(bond,4);
        %set the first neighbor for this particular vertex. this is not
        %particularly efficient: this will probably get reset a few times.
        %However, it is consistent.
        vdat(fbvert).nvert(1)=nbvert;
        vdat(fbvert).ncell(1)=ocell;
        vdat(fbvert).ncell(3)=i-1;
        
        %this section finds the corresponding bond number in the other cell
        %to help find neighboring vertices and cells in our vert structure.
        onbneighb=length(lat.cells{ocell+1});
        numflag=0;
        for k=1:onbneighb
            obond=lat.cells{ocell+1}(k);
            ofbvert=lat.bonds(obond,1);
            onbvert=lat.bonds(obond,2);
            if ofbvert==nbvert && onbvert==fbvert
                numflag=k;
            end
        end
        if numflag==0
            fprintf('ofbvert = %d, onbvert = %d, ',ofbvert,onbvert);
            fprintf('nbvert = %d, fbvert = %d\n',nbvert,fbvert);
            error('not good! logical error for numflag!\n');
        elseif numflag~=0
            bo_index=mod(numflag,onbneighb)+1;
        end
        b_index=mod(j-2,nbneighb)+1;
        %b_index and bo_index are the bond indices (in the corresponding
        %cells) that refer to bonds that connect to the vertex of interest
        %(in this case fbvert).  therefore the global bond indices should
        %be lat.cells{i+1}(b_index) for the first neighbor, while the
        %other should be lat.cells{ocell}(bo_index)
        bo_num=lat.cells{ocell+1}(bo_index);
        b_num=lat.cells{i}(b_index);
        
        %after the relevant bond indices are obtained, we can then use the
        %orientation of the bonds to determine the remaining neighboring
        %vertices.  recall that 1->2 in the bond index goes clockwise
        %around a cell according to their data organization. this part
        %gives the neighboring vertices. (draw a picture!).  note that i am
        %organizing the neighbors in counterclockwise order, and indexing
        %vertices <-> cells in the same fashion as the Hamiltonian code.
        vdat(fbvert).nvert(2)=lat.bonds(bo_num,2);
        vdat(fbvert).nvert(3)=lat.bonds(b_num,1);
        %this part gives the remaining neighboring cell, along with a
        %sanity check just to make sure things worked.
        if lat.bonds(bo_num,3)==ocell
            pcell1=lat.bonds(bo_num,4);
        elseif lat.bonds(bo_num,4)==ocell
            pcell1=lat.bonds(bo_num,3);
        else
            fprintf('this is bo_num=%d, with ocell=%d\n',bo_num,ocell);
            error('oops! somehow this bond does not border our ocell cell!\n');
        end
        if lat.bonds(b_num,3)==i-1
            pcell2=lat.bonds(b_num,4);
        elseif lat.bonds(b_num,4)==i-1
            pcell2=lat.bonds(b_num,3);
        else
            error('oops! somehow this bond does not border our i cell!\n');
        end
        if pcell1~=pcell2
            fprintf('pcell1 = %d, pcell2 = %d\n',pcell1,pcell2);
            fprintf('current cell = %d, other cell = %d\n',i,ocell);
            error('ouch. they do not agree on this neighboring cell!');
            %one thing to look out for potentially are more than 3-fold
            %vertices.  that would cause this to fail.  utterly. hopefully
            %their minimization program doesn't allow for that and this is
            %hopefully valid data for their minimization program.
        else
            %if all else has worked out, then here we declare our final
            %neighboring cell to this vertex and we are done!
            vdat(fbvert).ncell(2)=pcell1;
        end
    end
    %note that almost all of these vertices will be done and then redone
    %with each neighboring cell.  However, this is at most a factor of 3 in
    %a very short amount of time and it is only run once per time slice. so
    %this should not be too huge a time-waster hopefully!
end
        

for i=1:numvert
    vdat(i).index=i-1;
    vdat(i).xcoord=lat.verts(i,1)/25.0;
    vdat(i).ycoord=lat.verts(i,2)/25.0;
    %this is to change the vertex indexing so it matches alberto's format
    for j=1:3
        vdat(i).nvert(j)=vdat(i).nvert(j)-1;
    end
end
%this is to index the vertices in the same fashion as alberto's format
for i=1:length(cdat)
    nbneighb=cdat(i).nbneighb;
    for j=1:nbneighb
        cdat(i).nvert(j)=cdat(i).nvert(j)-1;
    end
end



function [newlat oldlat newcell0] = pare_data(lat,t)
% this pares the number bonds so that it is invertible by our inversion
% scheme. the t variable is just to help us identify what time this is
% happening at.  debugging variable.

numbond=length(lat.bonds);
bbonds=lat.cells{1}; %border bonds
nbbonds=length(bbonds);

oldlat=lat;
newlat=lat;

numnewbonds=length(newlat.bonds);
%initial procedure: start with first bond in list, and go "clockwise" to
%check if the neighboring cell also belongs to the same cell.  iterate
%until we find a different cell.  continue all the way around. note that
%"clockwise" for cell0 seem like it's c-clockwise, due to cell0 being
%the outside

%this is the list of cells bordering cell0 in accordance with the cell0
%index.  in other words, bclist(i) = {the cell that borders lat.cells{1}(i)
%that is not cell0}.  here k counts how many different cells we have on the
%border. and cell0list(k) tells us where the k'th cell change occurs
%(returns index i from bclist).
bclist(1:nbbonds)=-1; k=0;
for i=1:nbbonds;
    if lat.bonds(bbonds(i),3)~=0
        error('somehow i mistakenly assumed something about the structure');
    else
        ccell=lat.bonds(bbonds(i),4);
    end
    bclist(i)=ccell;
    % notice this next if statement only assigns values IF the neighboring
    % cell changes from the i-1 to the ith bond.
    if i>1 && bclist(i)~=bclist(i-1)
        k=k+1;
        cell0list(k)=i;
        shortbclist(k)=ccell; %short list of bordering cells
    end
end
if bclist(nbbonds)~=bclist(1)
    k=k+1;
    cell0list(k)=1;
    shortbclist(k)=lat.bonds(bbonds(1),4);
end
nbcellbound=k; %number of bounding cells to the system.

%now we can begin to create new bonds. note that we have to do this with
%cell0 as well as the cells that border cell0.  due to the way that their
%bond data structure works, it is unnecessary to delete the old bonds.

newcell0(1:nbcellbound)=-1; %this is the new cell0 (lat.cells{1}) array.
delvert(1:(nbbonds-nbcellbound))=-1; %this is the array of vertices that are deleted
delvert_counts=0;
for k=1:nbcellbound
    curr_bdcell_idx=shortbclist(k)+1;
    nbneighb=length(lat.cells{curr_bdcell_idx});
    bdcount=0;
    clear bd_idx; %crucial initialization for bd_idx thru loop
    for j=1:nbneighb
        % check if cell neighboring j'th bond is cell0
        curr_bond=lat.cells{curr_bdcell_idx}(j);
        if lat.bonds(curr_bond,4)==0
            bdcount=bdcount+1; %counts how many boundary bonds there are
            bd_idx(bdcount)=j; %list of boundary bonds
        end
    end
    if bdcount==0 %error checking
        error('supposedly no boundary bonds?  bad!\n');
    end
    bondfirst=lat.cells{curr_bdcell_idx}(1);
    bondlast=lat.cells{curr_bdcell_idx}(nbneighb);
    clear newcell
    newcell(1:nbneighb-bdcount+1)=0;
    if bdcount==1
        %basically nothing has to happen here. we leave it as is. it
        %already only borders cell0 once.
        for j=1:length(lat.cells{1})
            %find which original bond it is that satisfies our requirement
            if lat.bonds(lat.cells{1}(j),4)==(curr_bdcell_idx-1)
                newcell0(k)=lat.cells{1}(j); %set it to the appropriate bond
            end
        end
    elseif bdcount>1 && (lat.bonds(bondfirst,4)~=0 || lat.bonds(bondlast,4)~=0)
        %create new bond in the obvious way
        numnewbonds=numnewbonds+1;
        newlat.bonds(numnewbonds,1)=lat.bonds(lat.cells{curr_bdcell_idx}(bd_idx(1)),1);
        newlat.bonds(numnewbonds,2)=lat.bonds(lat.cells{curr_bdcell_idx}(bd_idx(bdcount)),2);
        newlat.bonds(numnewbonds,3)=curr_bdcell_idx-1;
        newlat.bonds(numnewbonds,4)=0;
        %create new cell structure
        for j=1:(nbneighb-bdcount+1)
            if j<bd_idx(1)
                newcell(j)=lat.cells{curr_bdcell_idx}(j);
            elseif j==bd_idx(1)
                newcell(j)=numnewbonds;
            elseif j>bd_idx(1)
                newcell(j)=lat.cells{curr_bdcell_idx}(j-1+bdcount);
            end
        end
        newlat.cells{curr_bdcell_idx}=newcell;
        
        %keeps track of which vertices need to be deleted.  to be taken
        %care of later
        for j=1:(bdcount-1)
            delvert_counts=delvert_counts+1;
            delvert(delvert_counts)=lat.bonds(lat.cells{curr_bdcell_idx}(bd_idx(j)),2);
        end
        
        %do the same stuff for cell0 now. create new bonds and put in the
        %new neighboring structure part
        numnewbonds=numnewbonds+1;
        newlat.bonds(numnewbonds,1)=newlat.bonds(numnewbonds-1,2);
        newlat.bonds(numnewbonds,2)=newlat.bonds(numnewbonds-1,1);
        newlat.bonds(numnewbonds,3)=0;
        newlat.bonds(numnewbonds,4)=curr_bdcell_idx-1;
        newcell0(k)=numnewbonds;
        
    elseif bdcount>1 && (lat.bonds(bondfirst,4)==0 && lat.bonds(bondlast,4)==0)
        %create new bond and cell structure in the not-so-obvious way.
        numnewbonds=numnewbonds+1;
        
        %this part searches for the "first" and "last" bonds (given our
        %orientation) that border cell0.
        prev=0;next=0;
        while 1==1
            if lat.bonds(lat.cells{curr_bdcell_idx}(nbneighb-prev),4)==0
                prev=prev+1;
            else
                break;
            end
        end
        while 1==1
            if lat.bonds(lat.cells{curr_bdcell_idx}(1+next),4)==0
                next=next+1;
            else
                break;
            end
        end
        if ((next+prev)~=bdcount)
            next
            prev
            bd_idx
            curr_bdcell_idx-1
            lat.cells{curr_bdcell_idx}
            error('bad counting for next and prev! check again!\n');
        end
        prev=prev-1;next=next-1; %this is just how i choose to rep. where the bonds are
        newlat.bonds(numnewbonds,1)=lat.bonds(lat.cells{curr_bdcell_idx}(nbneighb-prev),1);
        newlat.bonds(numnewbonds,2)=lat.bonds(lat.cells{curr_bdcell_idx}(1+next),2);
        newlat.bonds(numnewbonds,3)=curr_bdcell_idx-1;
        newlat.bonds(numnewbonds,4)=0;
        
        %create new cell structure
        newcell(1)=numnewbonds;
        for j=2:(nbneighb-bdcount+1)
            newcell(j)=lat.cells{curr_bdcell_idx}(j+next);
        end
        newlat.cells{curr_bdcell_idx}=newcell;
        
        %keeps track of which vertices need to be deleted.  to be taken
        %care of later
        for j=1:(bdcount-1)
            delvert_counts=delvert_counts+1;
            delvert(delvert_counts)=lat.bonds(lat.cells{curr_bdcell_idx}(mod(nbneighb-prev-2+j,nbneighb)+1),2);
        end
        
        %create for cell0 now. should be same as above part
        numnewbonds=numnewbonds+1;
        newlat.bonds(numnewbonds,1)=newlat.bonds(numnewbonds-1,2);
        newlat.bonds(numnewbonds,2)=newlat.bonds(numnewbonds-1,1);
        newlat.bonds(numnewbonds,3)=0;
        newlat.bonds(numnewbonds,4)=curr_bdcell_idx-1;
        newcell0(k)=numnewbonds;
    else
        error('logicall fallacy in bd_idx stuff!\n');
    end
end
newlat.cells{1}=newcell0;

if delvert_counts ~= (nbbonds-nbcellbound)
    delvert
    nbbonds
    nbcellbound
    shortbclist
    bclist
    t
    error('counting error!  need to check!\n');
end
delvert=sort(delvert);

%this makes the new vertex structure with the pared down number of vertices
count=0;newvert(1:(length(lat.verts)-length(delvert)),1:3)=-1;
for i=1:(length(lat.verts)-length(delvert))
    if count < length(delvert)
        while (i+count)==delvert(count+1)
            count=count+1;
            if count == length(delvert)
                break;
            end
        end
    end
    newvert(i,1)=lat.verts(i+count,1);
    newvert(i,2)=lat.verts(i+count,2);
    newvert(i,3)=lat.verts(i+count,3);
end
if count~=length(delvert)
    fprintf('count = %d, length(delvert) = %d\n',count,length(delvert));
    %error('error in counting!');
end
newlat.verts=newvert;

for i=1:length(newlat.bonds)
    v1=newlat.bonds(i,1);v2=newlat.bonds(i,2);
    count1=0;count2=0;
    while v1>delvert(count1+1)
        count1=count1+1;
        if count1==length(delvert)
            break;
        end
    end
    while v2>delvert(count2+1)
        count2=count2+1;
        if count2==length(delvert)
            break;
        end
    end
    newlat.bonds(i,1)=v1-count1;
    newlat.bonds(i,2)=v2-count2;
end

%everything up to here works for sure.  after this, we need to pare the
%vertex structure again, to get rid of vertices that we don't use. vertlist
%is the list of vertices that we do use. vcount counts the size of vertlist
numcell=length(newlat.cells);
clear vertlist;clear newvert;

for i=1:length(newlat.cells{1})
    vertlist(i)=newlat.bonds(newlat.cells{1}(i),1);
end


%this stuff expands vertlist
vcount=length(vertlist);
for i=2:numcell
    nbneighb=length(newlat.cells{i});
    for j=1:nbneighb
        vert=newlat.bonds(newlat.cells{i}(j),1);
        flag=0;
        for k=1:vcount
            if vert==vertlist(k)
                flag=1;
                break;
            end
        end
        if flag==0
            vcount=vcount+1;
            vertlist(vcount)=vert;
        end
    end
end
vertlist=sort(vertlist);

if vcount~=length(vertlist)
    error('bad counting!');
end
newvert(1:vcount,1:3)=0;
newbond=newlat.bonds;
for i=1:length(vertlist)
    newvert(i,:)=newlat.verts(vertlist(i),:);
    for j=1:length(newlat.bonds)
        if newlat.bonds(j,1)==vertlist(i)
            newbond(j,1)=i;
        end
        if newlat.bonds(j,2)==vertlist(i)
            newbond(j,2)=i;
        end
    end
end
newlat.bonds=newbond;newlat.verts=newvert;


%draws the new cell0. for testing purposes.
nbneighb=length(newlat.cells{1});
for i=1:nbneighb
    vert1x=newlat.verts(newlat.bonds(newlat.cells{1}(i),1),1);
    vert1y=newlat.verts(newlat.bonds(newlat.cells{1}(i),1),2);
    vert2x=newlat.verts(newlat.bonds(newlat.cells{1}(i),2),1);
    vert2y=newlat.verts(newlat.bonds(newlat.cells{1}(i),2),2);
    vx=[vert1x,vert2x]';vy=[vert1y,vert2y]';
    line(vx,vy,'Color',[0 0 0]);
end
