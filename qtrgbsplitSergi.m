% % The method uses three thresholds for red,green and blue blocks.
% % For any 3 channel block if the maximum value of the block elements minus the
% % minimum value of the block elements is greater than its respective threshold
% % then the 3 channel image is splitted into 4 blocks. 

% % Thresholds are specified as a value between 0 and 1, even if I is of class uint8 or
% % uint16. If I is uint8, the threshold value you supply is multiplied by
% % 255 to determine the actual threshold to use; if I is uint16, the
% % threshold value you supply is multiplied by 65535. If I is double (HSV
% % i.e.), no further change.

% % Input: image I (true color) image matrix of size M-by-N-by-3
% %        minimum block size desired on split process. 
% % Output: if image is split the method returns 
% %         size and starting coordinate indices of each block. If the
% %         image is not divided then it return -1 value for each output argument

function [b1r,b1c,b2r,b2c,b3r,b3c,b4r,b4c,szr,szc]=qtrgbsplitSergi(I,minBlockSize,varargin)

% % Default Values 
thvec=[0 0 0];                %lossless encoding
defaultValues = {thvec};      %threshold vector
% % % Assign Values
nonemptyIdx = ~cellfun('isempty',varargin);
defaultValues(nonemptyIdx) = varargin(nonemptyIdx);
[thvec] = deal(defaultValues{:});
% % %------------------------------------------------------
% % Handling Threshold vector for various
if(isempty(thvec))
    error('threshold vector must have atleast one value');
elseif(length(thvec)==1)
    ThR=thvec(1);
    ThG=ThR;
    ThB=ThG;
 elseif(length(thvec)==2)
    ThR=thvec(1);
    ThG=thvec(2);
    ThB=ThG;
 elseif(length(thvec)>=3)
    ThR=thvec(1);
    ThG=thvec(2);
    ThB=thvec(3);
end

if(ThR<0 || ThR>1 || ThG<0 || ThG>1 || ThB<0|| ThB>1)
    error('threshold value(s) must be between 0 and 1');
end

str=class(I);
if(strcmp(str,'uint8'))
    ThR=ThR*255;
    ThG=ThG*255;
    ThB=ThB*255;
elseif(strcmp(str,'uint16'))
    ThR=ThR*65535;
    ThG=ThG*65535;
    ThB=ThB*65535;
elseif(strcmp(str,'double'))
else    
    error('Input Image must be of type uint8 or uint16');
end
% % %------------------------------------------------------
% % % ThR,ThG,ThB

[ir ic id]=size(I);
if(id~=3)
    error('Input Matrix (Image) must be of size M-by-N-by-3 (true color)');    
end

b1r=-1;
b1c=-1;
b2r=-1;
b2c=-1;
b3r=-1;
b3c=-1;
b4r=-1;
b4c=-1;
szr=-1;
szc=-1;

if(ir==minBlockSize && ic==minBlockSize)  % block of size=1 are not split
    return
end

I=double(I); %convert to double for processing
stdLayer1 = std2(I(:,:,1));
stdLayer2 = std2(I(:,:,2));
stdLayer3 = std2(I(:,:,3));
splitCond = stdLayer1 > ThR | stdLayer2 > ThG | stdLayer3 > ThB;
if splitCond
   [b1r,b1c,b2r,b2c,b3r,b3c,b4r,b4c,szr,szc]=splitinto4(ir,ic);
end


% % % --------------------------------
% % % Author: Dr. Murtaza Khan
% % % Email : drkhanmurtaza@gmail.com
% % % Modification: Sergi Salgueiro Genís
% % % --------------------------------
