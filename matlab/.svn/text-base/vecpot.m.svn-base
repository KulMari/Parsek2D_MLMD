function ay=vecpot(tx,tz,tbx,tbz)

x=tx';z=tz';bx=tbx';bz=tbz';

[nx,nz]=size(x);

nzmezzo=ceil(nz/2);

ay=zeros(size(bx));
for ind=2:nx;
	ay(ind,nzmezzo)=ay(ind-1,nzmezzo)+trapz(x(ind-1:ind,nzmezzo),bz(ind-1:ind,nzmezzo));
end
for ind=nzmezzo+1:nz
%   ind;
   ay(:,ind)=ay(:,ind-1)-trapz(z(2,ind-1:ind),bx(:,ind-1:ind)')';
end
for ind=nzmezzo-1:-1:1
   ay(:,ind)=ay(:,ind+1)+trapz(z(2,ind:ind+1),bx(:,ind:ind+1)')';
end
ay=ay';