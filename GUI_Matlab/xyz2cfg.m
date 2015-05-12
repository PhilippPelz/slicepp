function xyz2cfg(filename)               % filename without extension

% written by Yaron Kauffmann - 01.11.08
% reads xyz files containing atomic positions in cartesian coordinates and converts to cfg format 
% readable by QSTEM (http://www.mf.mpg.de/en/organisation/hsm/koch/stem/index.html).
% the xyz file format should be as follows:
% line 1 : the number of atoms (natom)
% line 2 : empty or one integer number (no text)
% line 3 ... : atom type and x y z positions in cartesian coordinates
% SORTCELL function was downloaded from http://www.mathworks.de/matlabcentral/fileexchange/13770.
% writeCFG function was provided by Christoph Koch (developer of QSTEM).


%Read a "xyz" file

  name = ['H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si','P ','S ','Cl',...
    'Ar','K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',...
    'Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',...
    'I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',...
    'Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',...
    'Fr','Ra','Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr'];

  fid = fopen(strcat(filename,'.xyz'),'r');
  natom_temp = fscanf(fid,'%i');
  natom=natom_temp(1)

  for k1=1:natom
  data(k1,:)=textscan(fid,'%2c %f32 %f32 %f32',1);
  end

  
  b=sortcell(data,1);       % sort the data according to atom type

  coords= cell2mat(b(:,2:4));
  
  Mm=diag(max(coords)-min(coords)) ;
  coords=coords-repmat(min(coords),natom,1);

  for k1=1:natom
  aType(k1)=(1+findstr(char(b(k1,1)),name))/2;
  end
  
 writeCFG(filename,Mm,natom,aType,coords)
  
  fclose(fid);    



  
  
  function Y = sortcell(X, DIM)
% SORTCELL    Sort a cell array in ascending order.
%
% Description: SORTCELL sorts the input cell array according to the
%   dimensions (columns) specified by the user.
%
% Usage: Y = sortcell(X, DIM)
%
% Input:
%	   X: the cell array to be sorted.
%  DIM: (optional) one or more column numbers. Simply an array of one or
%       more column numbers.  The first number is the primary column on
%       which to sort. Extra column numbers may be supplied if secondary
%       sorting is required. The defuault value is 1, if no dimension
%       array is supplied.
%
% Output:
%     Y: the sorted cell array.
%
% Example:    Y = sortcell(X, [3 2])
%
% Note that this function has only been tested on mixed cell arrays
% containing character strings and numeric values.
% 
% Documentation Date: Feb.01,2007 13:36:04
% 
% Tags:
% {TAG} {TAG} {TAG}
% 
% 

%   Copyright 2007  Jeff Jackson (Ocean Sciences, DFO Canada)
%   Creation Date: Jan. 24, 2007
%   Last Updated:  Jan. 25, 2007

% Check to see if no input arguments were supplied.  If this is the case,
% stop execution and output an error message to the user.
if nargin == 0
	error('No input arguments were supplied.  At least one is expected.');
% Check to see if the only input argument is a cell array.  If it isn't
% then stop execution and output an error message to the user. Also set the
% DIM value since it was not supplied.
elseif nargin == 1
	if ~iscell(X)
		error('Input argument is not a cell array.  A cell array is expected.');
	end
	DIM = 1;
% Check to see if the first input argument is a cell array and the second
% one is numeric.  If either check fails then stop execution and output an
% error message to the user.
elseif nargin == 2
	if ~iscell(X)
		error('The first input argument is not a cell array.  A cell array is expected.');
	end
	if ~isnumeric(DIM)
		error('The second input argument is not numeric.  At least one numeric value is expected.');
	end
% Check to see if too many arguments were input.  If there were then exit
% the function issuing a error message to the user.
elseif nargin > 2
	error('Too many input arguments supplied.  Only two are allowed.');
end

% Now find out if the cell array is being sorted on more than one column.
% If it is then use recursion to call the sortcell function again to sort
% the less important columns first. Repeat calls to sortcell until only one
% column is left to be sorted. Then return the sorted cell array to the
% calling function to continue with the higher priority sorting.
ndim = length(DIM);
if ndim > 1
	col = DIM(2:end);
	X = sortcell(X, col);
end

% Get the dimensions of the input cell array.
[nrows, ncols] = size(X);

% Retrieve the primary dimension (column) to be sorted.
col = DIM(1);

% Place the cells for this column in variable 'B'.
B = X(:,col);

% Check each cell in cell array 'B' to see if it contains either a
% character string or numeric value. If it is a character string it returns
% a '1' in the same location of boolean array 'a'; a '0' otherwise. If it
% is a numeric value it returns a '1' in the boolean array 'b'; a '0'
% otherwise.
a = cellfun('isclass', B, 'char');
suma = sum(a);
b = cellfun('isclass', B, 'double');
sumb = sum(b);

% Check to see if cell array 'B' contained only character string.
% If cell array B contains character strings then do nothing because
% no further content handling is required.
if suma == nrows
% Check to see if cell array 'B' contained only numeric values.
elseif sumb == nrows
  % If the cells in cell array 'B' contain numeric values retrieve the cell
  % contents and change 'B' to a numeric array.
  B = [B{:}];
else
	error('This column is mixed so sorting cannot be completed.');
end

% Sort the current array and return the new index.
[ix,ix] = sort(B);

% Using the index from the sorted array, update the input cell array and
% return it.
Y = X(ix,:);


% function writeCFG(fileName,Mm,N,aType,coords,shift,mode)
%
% fileName:   name of cfg file
% Mm:         metric matrix (for orthogonal unit cells or super structures
%             Mm = diag([ax by cz]), where ax, by, cz are the size of your
%             unit cell / super structure in A.
% aType:      column vector of Z numbers
% coords:     N x 3 matrix of atomic positions [x1 y1 z1; x2 y2 z2; ...]
% shift:      can be used to define an offset to atomic positions (default
%             is [0 0 0]
% mode:       0 = atom positions in coords are fractional coordinates
%             1 = atom positions in coords are cartesian coordinates
function writeCFG(filename,Mm,N,aType,coords,shift,mode)

    mode=0;
    
if nargin < 6
    shift = [0 0 0];
end
if nargin < 7
    mode = 1;
end

Mminv = inv(Mm);
atomType = 0;
name = ['H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si','P ','S ','Cl',...
    'Ar','K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',...
    'Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',...
    'I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',...
    'Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',...
    'Fr','Ra','Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr'];
mass = 2*[1:length(name)];

fid = fopen(strcat(filename,'.cfg'),'w');

fprintf(fid,'Number of particles = %d\n',N);
fprintf(fid,'A = 1.0 Angstrom (basic length-scale)\n');
for ix=1:3
    for iy=1:3
        fprintf(fid,'H0(%d,%d) = %g A\n',ix,iy,Mm(ix,iy));
    end
end
fprintf(fid,'.NO_VELOCITY.\n');
fprintf(fid,'entry_count = 3\n');

for ix=1:N
    at = aType(ix);
    if at ~= atomType
        fprintf(fid,'%g\n%s\n',mass(at),name(2*at-1:2*at));  % specify type and mass of first atom (all others are the same)
        atomType = at;
    end
    if mode == 1
%        fprintf(fid,'%g %g %g\n',mod(coords(ix,:)*Mminv+shift,1));
        fprintf(fid,'%g %g %g\n',coords(ix,:)*Mminv+shift);
    else
        fprintf(fid,'%g %g %g\n',coords(ix,:));        
    end    
end
fclose(fid);