function M=select_sublat(L,fig)
%selects a sublattice from a larger bwlabel's image L

global Linternal
global cellvec
global done

figure(fig)
imshow(L);hold on;
impixelinfo;
Linternal=L;
cellvec=[];
set(gcf,'WindowButtonDownFcn',@wbd);
set(gcf,'WindowButtonUpFcn',@wbu);
set(gcf,'KeyPressFcn',@kpf);

done=0;
if(done==1)
    M=update_data(L,cellvec);
end


function wbd(obj,event)

global x0
global y0
if strcmp(get(gcf,'SelectionType'),'normal')
    pointer=get(gca,'CurrentPoint');
    x0=pointer(1,1);y0=pointer(1,2);
end

function wbu(obj,event)

global x0
global y0
global Linternal
global Lsub
global boxo

pointer=get(gca,'CurrentPoint');
xf=pointer(1,1);yf=pointer(1,2);
line([x0 xf xf x0 x0],[y0 y0 yf yf y0]);
xmin=min([x0,xf]);ymin=min([y0,yf]);xmax=max([x0,xf]);ymax=max([y0,yf]);
Lsub=Linternal(ymin:ymax,xmin:xmax);
%figure(11);imshow(Lsub);

function kpf(src,evnt)

global Lsub
global Linternal
global cellvec
global done

L=Linternal;

if evnt.Character=='A'
    %add stuff to the vector
    cellvecadd=unique(Lsub(Lsub>0));
    cellvec=unique([cellvec;cellvecadd]);
elseif evnt.Character=='C'
    %cancel adding this set
    boxo=findobj('Type','line');delete(boxo);
elseif evnt.Character=='D'
    %delete elements from vector
    Q=sprintf('Which cell(s) do you want to delete?  ');
    tmp=input(Q);
    cellvec=intersect(cellvec,setxor(cellvec,tmp));
    boxo=findobj('Type','line');delete(boxo);
elseif evnt.Character=='P'
    %print the vector
    cellvec
elseif evnt.Character=='U' %update
    zmod=update_data(L,cellvec);
elseif evnt.Character=='O' %output
    done=1;
elseif evnt.Character=='S' %save to .mat file
    [zmod,NewM1]=update_data(L,cellvec);
    save 'reducedlatticedata1.mat' NewM1;
    fprintf('saved!\n');
end

function [zmodfin,M]=update_data(L,cellvec)

fprintf('updating');
z=find(L==0);
fprintf('.');
sz=size(L);modfac=sz(1);
zmod=z(z>modfac & z<(sz(1)*sz(2)-modfac) & mod(z,modfac)>1);
hold on;
smallx=sz(2);smally=sz(1);bigx=0;bigy=0;count=0;
for i=1:length(zmod)
    if ismember(i,floor(length(zmod)/5*[1 2 3 4 5]))
        fprintf('.');
    end
    test=[L(zmod(i)-modfac-1),L(zmod(i)-modfac),L(zmod(i)-modfac+1),L(zmod(i)-1),L(zmod(i)),L(zmod(i)+1),L(zmod(i)+modfac-1),L(zmod(i)+modfac),L(zmod(i)+modfac+1)];
    if sum(ismember(test,cellvec))>0
        count=count+1;ZMODFIN(count)=zmod(i);
        currx=ceil(zmod(i)/modfac);curry=mod(zmod(i)-1,modfac)+1;
        if (currx+1>bigx)
            bigx=ceil(zmod(i)/modfac)+1;
        end
        if(currx-1<smallx)
            smallx=ceil(zmod(i)/modfac)-1;
        end
        if(curry+1>bigy)
            bigy=mod(zmod(i)-1,modfac)+2;
        end
        if(curry-1<smally)
            smally=mod(zmod(i)-1,modfac);
        end
        plot(currx,curry,'r.');
    end
end
zmodfin(1:length(ZMODFIN))=0;
zmodfin(:)=ZMODFIN(:);
imtmp(1:sz(1),1:sz(2))=0;
imtmp(zmodfin)=255;
figure(2);imshow(imtmp);
Lnew=watershed(imtmp);
M=Lnew(smally:bigy,smallx:bigx);
fprintf('done!\n');
