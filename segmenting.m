function [L,d1,E,Ux,Uy,I,G,J] = segmenting(E,Ux,Uy,I,G,J,l,B,l_k)


% all vectors must be sorted, starting element belongs to the most inner tube
%E, U, I, G, J   stifness, curvature, inertia, torsion constant, and second moment of inertia vectors for each tube
%l vector of tube length
%B  vector of tube movments with respect to template position, i.e., s=0 (always negative)
%l_k vecot oftube's curved part length

k=length(l); 

d1= l+B; % position of tip of the tubes
d2=d1-l_k; % position of the point where tube bending starts
points=[0 B d2 d1];
[L, index]=sort(points);
L = 1e-5*floor(1e5*diff(L));  % length of each segment 
%(used floor because diff command doesn't give absolute zero sometimes)

for i=1:k-1
if B(i)>B(i+1)
    sprintf('inner tube is clashing into outer tubes')
    E=zeros(k,length(L));
    I=E; G=E; J=E; Ux=E; Uy=E;
    return
end
end

EE=zeros(k,length(L));
II=EE; GG=EE; JJ=EE; UUx=EE; UUy=EE;

for i=1:k
    
a=find(index==i+1);   % find where tube begins
b=find(index==1*k+i+1); % find where tube curve starts
c=find(index==2*k+i+1); % find where tube ends

if L(a)==0; a=a+1;  end
if L(b)==0; b=b+1;  end
if c<=length(L)
    if L(c)==0; c=c+1; end
end
    
EE(i,a:c-1)=E(i);
II(i,a:c-1)=I(i);
GG(i,a:c-1)=G(i);
JJ(i,a:c-1)=J(i);
UUx(i,b:c-1)=Ux(i);
UUy(i,b:c-1)=Uy(i);
end

l=L(~(L==0));  % get rid of zero lengthes
E=zeros(k,length(l)); I=E; G=E; J=E; Ux=E; Uy=E;
 for i=1:k
    E(i,:)=EE(i,~(L==0)); I(i,:)=II(i,~(L==0)); G(i,:)=GG(i,~(L==0));
    J(i,:)=JJ(i,~(L==0)); Ux(i,:)=UUx(i,~(L==0)); Uy(i,:)=UUy(i,~(L==0));
 end
L=L(~(L==0));

end
