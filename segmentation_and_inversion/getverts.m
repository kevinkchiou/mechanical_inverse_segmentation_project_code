function [v c b vold]=getverts(L,mode)
%mode=0 implies we're taking the whole lattice.  use mode=1 for L's created
%with select_sublat() function

v=prelimverts(L,mode); %vertex detection in bwlabel'd image
[c v]=prelimcdat(v); %creation of vdat and cdat structures

[v c b]=prelimbdat(v,c,L); %creation of bond structures

%by this point all data structures should be created and nominally correct.
%this modifies bond and vertex indexing to agree with simulation indexing
[v c b]=indexmodification(v,c,b);
vold=v;

%this rescales lengths to something more appropriate
[v c b]=lengthrescale(v,c,b,20.0);


%%
function [v c b]=prelimbdat(v,c,L)
%creates bond structure

b(length(v)+length(c)-2).index=[];
nv=length(v);nb=0;
for i=1:nv
    for j=1:3
        k=v(i).nvert(j);
        if(k>i)
            nb=nb+1;
            b(nb).index=nb-1;
            b(nb).nvert=[i,k];
            b(nb).ncell(1)=v(i).ncell(j);b(nb).ncell(2)=v(i).ncell(mod(j,3)+1);
            b(nb).ncell(3)=v(i).ncell(mod(j+1,3)+1);
            l=find(v(k).nvert==i);
            if length(l)~=1
                fprintf('oops! we have i=%d,j=%d,k=%d, and l is wrong size.\n',i,j,k);
            end
            b(nb).ncell(4)=v(k).ncell(mod(l,3)+1);
            if v(k).ncell(l)~=b(nb).ncell(3) || v(k).ncell(mod(l+1,3)+1)~=b(nb).ncell(1)
                fprintf('seems to be some error in orientation or something!\n');
                fprintf('occurs in i = %d, and k = %d. check!\n',i,k);
            end
            x0=v(i).xcoord;y0=v(i).ycoord;
            xf=v(k).xcoord;yf=v(k).ycoord;
            %fprintf('working on bond %d...\n',nb);
            if(b(nb).ncell(1)~=0 && b(nb).ncell(3)~=0)
                [b(nb).pix,prob]=getbondpixels([x0;y0],[xf;yf],L);
            else
                b(nb).pix=[];prob=0;
            end
            if prob==1
                fprintf('we have i=%d,j=%d,k=%d,nb=%d\n',i,j,k,nb);
            end
        end
        %fprintf('going...');
    end
end
bold=b;

%created for useful vector operations
tempb(1:2,1:length(b))=0;
for i=1:length(b)
    tempb(1:2,i)=b(i).nvert(1:2)';
end

for i=1:nv
    cv=i;
    for j=1:3
        nv=v(i).nvert(j);
        if(cv>nv)
            test=tempb-repmat([nv;cv],1,length(tempb));
            val=find(test(1,:).^2 + test(2,:).^2 ==0);
            v(i).nbond(j)=val;
        elseif(nv>cv)
            test=tempb-repmat([cv;nv],1,length(tempb));
            val=find(test(1,:).^2 + test(2,:).^2 ==0);
            v(i).nbond(j)=val;
        end
    end
end
%hopefully will create less fragmented structure.
%{
bnew(nb).index=[];bnew(nb).nvert=[];bnew(nb).ncell=[];bnew(nb).nbond=[];
bnew(nb).length=[];bnew(nb).tension=[];bnew(nb).pix=[];
for i=1:nb
    bnew(i)=b(i);
end
clear b;
b=bnew;
%}

nb=length(b);
for i=1:nb
    fc=b(i).ncell(1)+1;
    rc=b(i).ncell(3)+1;
    fidx=find(c(fc).nvert==b(i).nvert(1));ridx=find(c(rc).nvert==b(i).nvert(2));
    c(fc).nbond(fidx)=i;c(fc).ncell(fidx)=rc-1;
    c(rc).nbond(ridx)=i;c(rc).ncell(ridx)=fc-1;
end
%%
function [c v vold]=prelimcdat(v)
%creates preliminary cdat structure with oriented vertices.

vc=v(:,3:5); %cell neighbors for each vertex.
vv(1:length(v),1:3)=0; %vertex neighbors.
numcell=max(max(vc));
c(numcell+1).index=numcell; %initializing structure
c(numcell+1).centx=[];c(numcell+1).centy=[];c(numcell+1).centz=0;
c(numcell+1).nvert=[];c(numcell+1).ncell=[];c(numcell+1).nbond=[];c(numcell+1).pressure=0;
c(numcell+1).nbneighb=0;
for i=1:numcell+1
    c(i).index=i-1;c(i).pressure=0;
    temp=mod(find(vc==c(i).index)-1,length(vc))+1;%index of verts neighboring cell i-1
    if(isempty(temp))
        fprintf('cell %d should be marked for deletion\n',i)
    else
        c(i).centx=mean(v(temp,1));c(i).centy=mean(v(temp,2));
        xpos=v(temp,1)-c(i).centx;ypos=v(temp,2)-c(i).centy;c(i).centz=0;
        [order,index]=sort(atan2(ypos,xpos));
        if i==1 %difference for cell0
            ordtemp=temp(index(end:-1:1)); %backwards from other cells.
            c(i).nvert=ordtemp';
        else
            ordtemp=temp(index);%put neighbors into structure element with orientation.
            c(i).nvert=ordtemp';
        end
        c(i).nbneighb=length(ordtemp);
        c(i).ncell=repmat(-1,1,length(ordtemp));
        c(i).nbond=repmat(-1,1,length(ordtemp));
    end
end
%need this for non-convex cell0
for i=1:c(1).nbneighb-1
    %cidx=neighb cell index "previous" to cell0 for this vertex
    cidx=vc(c(1).nvert(i),mod(find(vc(c(1).nvert(i),:)==0)+1,3)+1);
    %this finds the "previous" vertex in that cell to the current vertex.
    %the net result is that this "previous" vertex is "next" for cell0.
    c(1).nvert(i+1)=c(cidx+1).nvert(mod(find(c(cidx+1).nvert==c(1).nvert(i))-2,length(c(cidx+1).nvert))+1);
end
%once they're in order determine vv topology

for i=1:length(v)
    vv(i,1)=c(v(i,3)+1).nvert(mod(find([c(v(i,3)+1).nvert]==i),length([c(v(i,3)+1).nvert]))+1);
    vv(i,2)=c(v(i,4)+1).nvert(mod(find([c(v(i,4)+1).nvert]==i),length([c(v(i,4)+1).nvert]))+1);
    vv(i,3)=c(v(i,5)+1).nvert(mod(find([c(v(i,5)+1).nvert]==i),length([c(v(i,5)+1).nvert]))+1);
end
v(:,6:8)=vv(:,1:3); %add it on to v.
vold(1:length(v),1:8)=0;
vold=v; %save unfragmented v.

nvert=length(v);
vdat(nvert).xcoord=0;vdat(nvert).ycoord=0;vdat(nvert).zcoord=0.0;
vdat(nvert).nvert=[-1 -1 -1];vdat(nvert).ncell=[-1 -1 -1];vdat(nvert).nbond=[-1 -1 -1];
vdat(nvert).xforce=0;vdat(nvert).yforce=0;vdat(nvert).zforce=0;
vdat(nvert).vdistx=[0 0 0];vdat(nvert).vdisty=[0 0 0];vdat(nvert).vdist=[0 0 0];
for i=1:nvert
    vdat(i).xcoord=v(i,1);vdat(i).ycoord=v(i,2);vdat(i).zcoord=0.0;
    vdat(i).nvert=v(i,6:8);vdat(i).ncell=v(i,3:5);
    vdat(i).index=i-1;
    vdat(i).vdistx=v(v(i,6:8),1)'-vdat(i).xcoord;
    vdat(i).vdisty=v(v(i,6:8),2)'-vdat(i).ycoord;
    vdat(i).vdist=sqrt([vdat(i).vdistx].^2 + [vdat(i).vdisty].^2);
    vdat(i).xforce=0;vdat(i).yforce=0;vdat(i).zforce=0;
end
clear v;
v=vdat;%output of new vertex structure

%%
function v=prelimverts(L,mode)
%takes a bwlabel'd image L and returns vertex locations v and the cells it
%neighbors as a vector. image processing step to identify relevant vertices
sz=size(L);
modfac=sz(1);

z=find(L==0); %vector of zeros in L
%these are elements of z that aren't along the edges of the picture. this
%is important since we want to use a 3x3 square to detect vertices.
%zmod = pixels where L==0 that are NOT along the border of image.
zmod=z(z>modfac & z<(sz(1)*sz(2)-modfac) & mod(z,modfac)>1);

j=0;
for i=1:length(zmod)
    %list of all pixels surrounding the one of interest
    test=[L(zmod(i)-modfac-1),L(zmod(i)-modfac),L(zmod(i)-modfac+1),L(zmod(i)-1),L(zmod(i)),L(zmod(i)+1),L(zmod(i)+modfac-1),L(zmod(i)+modfac),L(zmod(i)+modfac+1)];
    testvec=unique(test);
    if(length(testvec)==4 && testvec(1)==0)
        j=j+1;
        v(j,1)=ceil(zmod(i)/modfac);v(j,2)=mod(zmod(i)-1,modfac)+1;
        v(j,3:5)=testvec(2:4)';
    end
end

j=0;
for i=1:size(v)
    if(v(i,3)>0 && v(i,4)>0 && v(i,5)>0)
        j=j+1;
        vnew(j,:)=v(i,:);
    end
end

s=regionprops(L);
cent=cat(1,s.Centroid);
v=vnew;vnew(:,3:5)=-1;
for i=1:size(v);
    if(v(i,3)==1)
        cent3=1.2*(v(i,1:2)-cent(v(i,3),:))+cent(v(i,3),:);
    else
        cent3=cent(v(i,3),:);
    end
    if(v(i,4)==1)
        cent4=1.2*(v(i,1:2)-cent(v(i,4),:))+cent(v(i,4),:);
    else
        cent4=cent(v(i,4),:);
    end
    if(v(i,5)==1)
        cent5=1.2*(v(i,1:2)-cent(v(i,5),:))+cent(v(i,5),:);
    else
        cent5=cent(v(i,5),:);
    end
    temp(1,:)=cent3-v(i,1:2);
    temp(2,:)=cent4-v(i,1:2);
    temp(3,:)=cent5-v(i,1:2);

    otemp(1)=atan2(temp(1,2),temp(1,1));
    otemp(2)=atan2(temp(2,2),temp(2,1));
    otemp(3)=atan2(temp(3,2),temp(3,1));
    [x,n]=sort(otemp);
    
    vnew(i,3)=v(i,n(1)+2);
    vnew(i,4)=v(i,n(2)+2);
    vnew(i,5)=v(i,n(3)+2);
end

if mode==0
    %pick up vertices along the edges
    zmod=z(mod(z,modfac)<2 | z<modfac+1 | z >(sz(1)*sz(2)-modfac));
    len=length(vnew);
    for i=1:length(zmod)
        vnew(len+i,1)=ceil(zmod(i)/modfac);vnew(len+i,2)=mod(zmod(i)-1,modfac)+1;
        if mod(zmod(i),modfac)==1
            vnew(len+i,3)=L(zmod(i)+modfac);vnew(len+i,4)=L(zmod(i)-modfac);
            vnew(len+i,5)=0;
        elseif mod(zmod(i),modfac)==0
            vnew(len+i,3)=L(zmod(i)-modfac);vnew(len+i,4)=L(zmod(i)+modfac);
            vnew(len+i,5)=0;
        elseif zmod(i)<modfac+1
           vnew(len+i,3)=L(zmod(i)-1);vnew(len+i,4)=L(zmod(i)+1);
           vnew(len+i,5)=0;
        elseif zmod(i)>sz(1)*sz(2)-modfac
            vnew(len+i,3)=L(zmod(i)+1);vnew(len+i,4)=L(zmod(i)-1);
            vnew(len+i,5)=0;
        end
    end
elseif mode==1
    j=0;
    vnewest=[];cellvec=2:1:max(max(L));
    for i=1:length(vnew)
        cells=vnew(i,3:5);
        membercelltest=ismember(cells,cellvec);
        if sum(~membercelltest)<2
            j=j+1;
            cells(~membercelltest)=0;
            cells(cells>0)=cells(cells>0)-1;
            vnew(i,3:5)=cells;
            vnewest(j,1:5)=vnew(i,1:5);
            %{
            if sum(vnewest(j,3:5)==0)>0
                vnewest(j,5:-1:3)=vnewest(j,3:5);
            end
            %}
        end
    end
    clear vnew;
    vnew=vnewest;
end

clear v;
l=size(vnew);
temp(l(1),l(2))=0;
temp(:,:)=vnew(:,:); %create a matrix that isn't so fragmented
v=temp;
%%
function [bp,chk]=getbondpixels(r1,r2,L)
%returns pixels involved in the bond between vertex positions r1 and r2

chk=0;
if size(r1)==[1 2]
    r1=r1';
end
if size(r2)==[1 2]
    r2=r2';
end

if (L(r1(2),r1(1))~=0 || L(r2(2),r2(1))~=0)
    error('we are not on a vertex!');
end
sz=size(L);edger1=0;edger2=0;
if(r1(1)==1 || r1(1)==sz(2))
    edger1=1;
end
if(r2(1)==1 || r2(1)==sz(2))
    edger2=1;
end
if(r1(2)==1 || r1(2)==sz(1))
    edger1=1;
end
if(r2(2)==1 || r2(2)==sz(1))
    edger2=1;
end

if edger1==1 && edger2==1
    bp=[];
    return;
end

v=r2-r1;

%pick direction of initial motion. dir is most obvious. diro is other
%choice. dirl is last choice.
testpos=[r1 r1 r1 r1]+[1 0 -1 0;0 1 0 -1];
dirvec=testpos(:,diag(L(testpos(2,:),testpos(1,:))==0));
if length(dirvec)>3
    %dirvec
    error('hmmm...');
end
dir1=dirvec(:,1)-r1;dir2=dirvec(:,2)-r1;dir3=dirvec(:,3)-r1;
dir1test=dir1'*v;dir2test=dir2'*v;dir3test=dir3'*v;
[srt sidx]=sort([dir1test dir2test dir3test]);
dirl=dirvec(:,sidx(1))-r1;diro=dirvec(:,sidx(2))-r1;dir=dirvec(:,sidx(3))-r1;
%{
if abs(v(1))>abs(v(2))
    dir=[sign(v(1));0];
    diro=[0;sign(v(2))];
    dirl=[0;-sign(v(2))];
    if L(r1(2)+dir(2),r1(1)+dir(1))~=0
        dir=[0;sign(v(2))];
        diro=[sign(v(1));0];
        dirl=[-sign(v(1));0];
    end
else
    dir=[0;sign(v(2))];
    diro=[sign(v(1));0];
    dirl=[-sign(v(1));0];
    if L(r1(2)+dir(2),r1(1)+dir(1))~=0
        dir=[sign(v(1));0];
        diro=[0;sign(v(2))];
        dirl=[0;-sign(v(2))];
    end
end
%}
testdir=[[1 0 -1 0];[0 1 0 -1];];diroflag=0;dirlflag=0;
ppos=r1;pos=r1+dir;npos=pos;step=0;bp(1:2,1)=ppos;bp(1:2,2)=pos;
while (step<(3*sqrt(v'*v)+1) && (npos(1)~=r2(1) || npos(2)~=r2(2)))
    testpos=[pos pos pos pos]+testdir;
    testpos=testpos(:,(testpos(1,:)>0 & testpos(2,:)>0)); %eliminates edges from testing
    testvec=diag(L(testpos(2,:),testpos(1,:)));
    npos=testpos(:,(testvec==0 & (((testpos(1,:)-ppos(1)).^2 + (testpos(2,:)-ppos(2)).^2)'~=0)));
    step=step+1; %counter
    %now reset ppos and pos for next iteration
    ppos=pos;pos=npos;
    %set the bond pixel
    if((sum(size(pos))>3 || step>(3*sqrt(v'*v)-2)) && diroflag==0) %reach a different vertex than anticipated! re-initialize!
        step=0;ppos=r1;pos=r1+diro;npos=pos;clear bp;bp(1:2,1)=ppos;bp(1:2,2)=pos;diroflag=1;
        %fprintf('doing diro in the (%d,%d) direction\n',diro(1),diro(2));
    elseif((sum(size(pos))>3 || step>(3*sqrt(v'*v)-2))&& diroflag==1 && dirlflag==0)
        step=0;ppos=r1;pos=r1+dirl;npos=pos;clear bp;bp(1:2,1)=ppos;bp(1:2,2)=pos;dirlflag=1;
        %fprintf('doing dirl in the (%d,%d) direction\n',dirl(1),dirl(2));
    else %not yet a vertex, have not reached limit: increment as usual.
        bp(:,step+2)=pos;
    end
end
if step>(3*sqrt(v'*v)-2)
    fprintf('we probably followed the wrong path. check out these particular vertices\n');
    chk=1;
end

%just to make sure final vector is fragmented as little as possible
bpfinal(1:2,1:step+2)=0; %initialize
bpfinal(:,:)=bp(:,:); %copy over term by term
clear bp; %clear the loop-constructed variable
bp=bpfinal';

%%
function [v c b]=indexmodification(vold,cold,bold)

v=vold;b=bold;c=cold;

for i=1:length(v)
    v(i).nvert=v(i).nvert-1;
    v(i).nbond=v(i).nbond-1;
end
for i=1:length(b)
    b(i).nvert=b(i).nvert-1;
end
for i=1:length(c)
    c(i).nvert=c(i).nvert-1;
    c(i).nbond=c(i).nbond-1;
end

cnew(1:length(c)-1)=c(2:length(c));
c=cnew;

function [v c b]=lengthrescale(vold,cold,bold,r)
%r is the scale factor

v=vold;b=bold;c=cold;
for i=1:length(v)
    v(i).xcoord=v(i).xcoord/r;
    v(i).ycoord=v(i).ycoord/r;
    v(i).vdistx=v(i).vdistx/r;
    v(i).vdisty=v(i).vdisty/r;
    v(i).vdist=v(i).vdist/r;
end
for i=1:length(b)
    if ~isempty(b(i).pix)
        b(i).pix=b(i).pix/r;
    end
end
for i=1:length(c)
    c(i).centx=c(i).centx/r;
    c(i).centy=c(i).centy/r;
end