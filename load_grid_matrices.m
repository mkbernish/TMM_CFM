%% load grids and matrices

load(matrix_path,'M3d','grid','MSKS',Amatrix); 
eval(['A = ' Amatrix ';']); eval(['clear ' Amatrix]);
m = size(A,1);
[j1, kprod] = min(abs((grid.zw+grid.dzt)-zprod))

% indices
iocn = find(M3d(:)==1); 
i3d = M3d*0+NaN; i3d(iocn) = 1:m;
isurf = find(M3d(:,:,1)==1);
[j1,ivsurf,j2] = intersect(iocn,isurf);
[j1,isub,j2] = intersect(iocn,find(M3d(:,:,1:kprod)==1)); % productivity subregion
[j1,idif] = setdiff(iocn,find(M3d(:,:,1:kprod)==1)); % no productivity subregion
[j1,icont] = intersect(iocn,find(grid.ZT3d<=100));

% volume and area
VOL = grid.DXT3d.*grid.DYT3d.*grid.DZT3d;
AREA = grid.DXT3d.*grid.DYT3d;
dv = VOL(iocn);
da = AREA(iocn);