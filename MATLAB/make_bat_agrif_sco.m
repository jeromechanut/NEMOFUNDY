% Program to ensure bathymetries matching 
% with AGRIF in NEMO s coordinate case (only, not s-z)
% Odd refinement factors only
% J. Chanut, September 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
path(path,'/homelocal-px/px-105/jchanut/MATLAB/LP_Bathymetry/Mfiles')
%
% Input bathymetries:
filebat_lev0_in='Bathymetry_BoF180_MITCH_NEW_range_0.5_7m_StJohnZZ15m_RedVol_sco.nc'; 		% Parent
%filebat_lev1_in='Bathymetry_SJAP100_NEW_range_0.5_7m_smooth.nc';
filebat_lev1_in='tmp.nc';		                        % Child
%filebat_lev1_in='/home/jchanut/BoF_DATA/Bathymetry_SJAP100_range_0.5_7m_spencer_smooth.nc';		                        % Child
%
% Output bathymetries:
filebat_lev0_out='Bathymetry_BoF180_MITCH_NEW_range_0.5_7m_StJohnZZ15m_RedVol_sco_agrif_new.nc';     % Parent
filebat_lev1_out='Bathymetry_SJAP100_NEW_range_0.5_7m_smooth_sco_agrif_new.nc';	% Child
%
fraf = 5;    % Refinement factor
rmax = 0.2;  % Maximum slope factor
nmatch = 3;  % Nb of matching points (in coarse grid points)
ntrans = 3;  % Nb of transition points (in coarse grid points)
%
% agrif Parameters as in input ascii file:
%
i0ag=374;
i1ag=462;
j0ag=191;
j1ag=362;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
i0ag_0=341;
i1ag_0=500;
j0ag_0=161;
j1ag_0=362;
%
% First and last interior coincident T-points on parent grid:
% Do not change
imin1=i0ag+2;
imax1=i1ag+1;
jmin1=j0ag+2;
jmax1=j1ag+1;
%
ntrans = ntrans*fraf; % Transition width in fine grid points
%
% 0. Read Bathymetries
%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%
%
%ncid = netcdf.open(filebat_lev0_in,'NOWRITE');
%varid = netcdf.inqVarID(ncid,'nav_lon');
%lon0 = netcdf.getVar(ncid,varid,'double');
%varid = netcdf.inqVarID(ncid,'nav_lat');
%lat0 = netcdf.getVar(ncid,varid,'double');
%varid = netcdf.inqVarID(ncid,'Bathymetry');
%bat0 = netcdf.getVar(ncid,varid,'double');
%netcdf.close(ncid);
nc = netcdf(filebat_lev0_in, 'nowrite');
lon0=nc{'nav_lon'}(:,:);
lat0=nc{'nav_lat'}(:,:);
bat0=nc{'Bathymetry'}(:,:);
nc = close(nc);
lon0=lon0';
lat0=lat0';
bat0=bat0';
[jpi0,jpj0]=size(bat0);
%
%ncid = netcdf.open(filebat_lev1_in,'NOWRITE');
%varid = netcdf.inqVarID(ncid,'nav_lon');
%lon1 = netcdf.getVar(ncid,varid,'double');
%varid = netcdf.inqVarID(ncid,'nav_lat');
%lat1 = netcdf.getVar(ncid,varid,'double');
%varid = netcdf.inqVarID(ncid,'Bathymetry');
%bat1 = netcdf.getVar(ncid,varid,'double');
%netcdf.close(ncid);
nc = netcdf(filebat_lev1_in, 'nowrite');
lon1=nc{'nav_lon'}(:,:);
lat1=nc{'nav_lat'}(:,:);
bat1=nc{'Bathymetry'}(:,:);
nc = close(nc);
lon1=lon1';
lat1=lat1';
bat1=bat1';
[jpi1,jpj1]=size(bat1);
%
% Crop child domain to match what is expected on AGRIF child grid domain:
% That's cooking here... 
i1=(i0ag-i0ag_0)*fraf+2;
% Line below because of 2 aditionnal lines in input child bathymetry
i2=jpi1-(i1ag_0-i1ag)*fraf-2;
%i2=jpi1-(i1ag_0-i1ag)*fraf;
%%%%%%%%%%%
j1=(j0ag-j0ag_0)*fraf+2;
j2=jpj1-(j1ag_0-j1ag)*fraf;
lon1=lon1(i1:i2,j1:j2);
lat1=lat1(i1:i2,j1:j2);
bat1=bat1(i1:i2,j1:j2);
[jpi1,jpj1]=size(bat1);
%
mask1=ones(size(bat1));
I=find(bat1==0);
mask1(I)=0.;
%
%
% Check first coincident point:
ifirst=3+(fraf-1)/2; % On child grid
'Child ', lon1(ifirst,ifirst), lat1(ifirst,ifirst)
'Parent', lon0(imin1,jmin1), lat0(imin1,jmin1)
%
%
bat2=bat1;
for ji=2:jpi1-1
    for jj=2:jpj1-1  
       tmp =  4.* mask1(ji,jj) + ...
              2.*(mask1(ji-1,jj  ) + mask1(ji+1,jj)    + ...
                  mask1(ji  ,jj-1) + mask1(ji,jj+1))   + ...
              1.*(mask1(ji-1,jj-1) + mask1(ji-1,jj+1)  + ...
                  mask1(ji+1,jj-1) + mask1(ji+1,jj+1));
       bat2(ji,jj) =  (4.* bat1(ji,jj) + ...
                       2.*(bat1(ji-1,jj  ) + bat1(ji+1,jj)    + ...
                           bat1(ji  ,jj-1) + bat1(ji,jj+1))   + ...
                       1.*(bat1(ji-1,jj-1) + bat1(ji-1,jj+1)  + ...
                           bat1(ji+1,jj-1) + bat1(ji+1,jj+1)) )/max(tmp,1.)*mask1(ji,jj);                 
    end
end 
bat1=GRID_PlusMinusScheme_rx0(mask1,bat2,rmax);
%
% Build weight array:
wgt1=bat1*0.;

% West:
%for j=1:jpj1
for j=1:171
   ist=2+nmatch*fraf;
   for l=1:ntrans
      wgt1(ist+l,j) = (ntrans-l)/ntrans ;
   end
end

% East:
%for j=1:jpj1
for j=1:400
   ist=jpi1-(2+nmatch*fraf)+1;
   for l=1:ntrans
      wgt1(ist-l,j) = max(wgt1(ist-l,j),(ntrans-l)/ntrans) ;
   end
end

% South:
for i=1:jpi1
   jst=2+nmatch*fraf;
   for l=1:ntrans
      wgt1(i,jst+l) = max(wgt1(i,jst+l),(ntrans-l)/ntrans) ;
   end
end

% North:
%for i=1:jpi1
%   jst=jpj1-(2+nmatch*fraf)+1;
%   for l=1:ntrans
%      wgt1(i,jst-l) = max(wgt1(i,jst-l),(ntrans-l)/ntrans) ;
%   end
%end

%for j=1:jpj1
%   wgt1(1:2+nmatch*fraf,j) = 1. ;
%   wgt1(jpi1-(2+nmatch*fraf)+1:jpi1,j) = 1. ;
%end
wgt1(1:2+nmatch*fraf,1:400) = 1. ;
wgt1(jpi1-(2+nmatch*fraf)+1:jpi1,1:400) = 1. ;
for i=1:jpi1
   wgt1(i,1:2+nmatch*fraf) = 1. ;
   wgt1(i,jpj1-(2+nmatch*fraf)+1:jpj1) = 1. ;
end
%
% Sort of "nearest neighbour" interpolation of parent grid bathymetry:
dx=(fraf-1)/2;
%bat00=bat1;
for i0=imin1-1:imax1+1
for j0=jmin1-1:jmax1+1
   i1 = 3+dx-fraf+(i0-(imin1-1))*fraf;
   j1 = 3+dx-fraf+(j0-(jmin1-1))*fraf;
   ii0=min(max(1,i1-dx),jpi1);
   ii1=min(max(1,i1+dx),jpi1);
   ji0=min(max(1,j1-dx),jpj1);
   ji1=micd 
n(max(1,j1+dx),jpj1);
   bat00(ii0:ii1,ji0:ji1) = bat0(i0,j0);
end
end
%
% Ensure bathymetry matching:
I=find(wgt1==1);
bat1(I)=bat00(I);
mask1=ones(size(bat1));
I=find(bat1==0);
mask1(I)=0.;
%
%
% 1. From highest resolution grid down to lowest, update bathymetry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
for i0=imin1:imax1
for j0=jmin1:jmax1
   i1 = 3+dx-fraf+(i0-(imin1-1))*fraf;
   j1 = 3+dx-fraf+(j0-(jmin1-1))*fraf;
   ii0=i1-dx;
   ii1=i1+dx;
   ji0=j1-dx;
   ji1=j1+dx;
   zbat=sum(sum(bat1(ii0:ii1,ji0:ji1)));
   ncount=sum(sum(mask1(ii0:ii1,ji0:ji1)));
   if (ncount<=0.5*fraf*fraf)
%   if (ncount==0)
      bat0(i0,j0)=0.;
   else
      bat0(i0,j0)=zbat/ncount;
%      bat0(i0,j0)=zbat/(fraf*fraf);
   end
end
end
%
mask0=ones(size(bat0));
I=find(bat0==0);
mask0(I)=0.;
%
% WRITE ouput
% Parent grid:
%nccreate(filebat_lev0_out,'nav_lon','dimensions',{'y' jpi0 'x' jpj0},'format','classic','Datatype','single')
%ncwrite(filebat_lev0_out,'nav_lon',lon0);
%nccreate(filebat_lev0_out,'nav_lat','dimensions',{'y' jpi0 'x' jpj0},'format','classic','Datatype','single')
%ncwrite(filebat_lev0_out,'nav_lat',lat0);
%nccreate(filebat_lev0_out,'Bathymetry','dimensions',{'y' jpi0 'x' jpj0},'format','classic','Datatype','single')
%ncwrite(filebat_lev0_out,'Bathymetry',bat0);
%
f = netcdf(filebat_lev0_out, 'clobber');
f.source = filebat_lev0_in;
f.version = '2.0';
f.Author = 'J. Chanut, Mercator Ocean';
f.Created = datestr(now);
f('x') = jpi0;
f('y') = jpj0;
f{'nav_lon'} = ncfloat({'y';'x'});
f{'nav_lon'}.long_name='Longitude';
f{'nav_lat'} = ncfloat({'y';'x'});
f{'nav_lat'}.long_name='Latitude';
%
%
%
f{'Bathymetry'} = ncfloat({'y';'x'});
f{'Bathymetry'}.long_name='Bathymetry';
f{'Bathymetry'}.units='m';
f{'Bathymetry'}.missing_value=ncfloat(0.);
%

f{'nav_lon'}(:,:)=lon0';
f{'nav_lat'}(:,:)=lat0';
f{'Bathymetry'}(:,:) = bat0';
%
close(f);
% Child grid:
%nccreate(filebat_lev1_out,'nav_lon','dimensions',{'y' jpi1 'x' jpj1},'format','classic','Datatype','single')
%ncwrite(filebat_lev1_out,'nav_lon',lon1);
%nccreate(filebat_lev1_out,'nav_lat','dimensions',{'y' jpi1 'x' jpj1},'format','classic','Datatype','single')
%ncwrite(filebat_lev1_out,'nav_lat',lat1);
%nccreate(filebat_lev1_out,'Bathymetry','dimensions',{'y' jpi1 'x' jpj1},'format','classic','Datatype','single')
%ncwrite(filebat_lev1_out,'Bathymetry',bat1);
%
f = netcdf(filebat_lev1_out, 'clobber');
f.source = filebat_lev1_in;
f.version = '2.0';
f.Author = 'J. Chanut, Mercator Ocean';
f.Created = datestr(now);
f('x') = jpi1;
f('y') = jpj1;
f{'nav_lon'} = ncfloat({'y';'x'});
f{'nav_lon'}.long_name='Longitude';
f{'nav_lat'} = ncfloat({'y';'x'});
f{'nav_lat'}.long_name='Latitude';
%
%
%
f{'Bathymetry'} = ncfloat({'y';'x'});
f{'Bathymetry'}.long_name='Bathymetry';
f{'Bathymetry'}.units='m';
f{'Bathymetry'}.missing_value=ncfloat(0.);
%

f{'nav_lon'}(:,:)=lon1';
f{'nav_lat'}(:,:)=lat1';
f{'Bathymetry'}(:,:) = bat1';
%
close(f);


