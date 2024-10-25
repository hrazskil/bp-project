function [eM2] = magFieldM(rObserved,dip,f,complAmpl)
%MAGFIELDM Summary of this function goes here
%inicialization
construct   = utilities.constants.giveConstants;
nDip    = size(dip.pos,1);
nObs    = size(rObserved,1);
eM2     = nan(nObs,nDip,3);     % Preallocate
MatrixRVec  = zeros(nObs,nDip,3);   % Preallocate
MR      = zeros(nObs,nDip);     % Preallocate
MdirR   = zeros(nObs,nDip,3);   % Preallocate
omega   = 2*pi*f;
k       = omega/construct.c0;


MatrixRVec  = reshape( repmat(rObserved,nDip,1) ,nObs, [],3) -...
    permute( reshape( repmat(dip.pos,nObs,1) ,nDip,[],3) ,[2 1 3]);
MR      = utilities.rowNormAccrossM(MatrixRVec);
MdirR   = MatrixRVec./MR;
MDDir   = permute( reshape( repmat(dip.dir,nObs,1) ,nDip,[],3) ,[2 1 3]);
Mp      = repmat(complAmpl',nObs,1,3) .* ( ...
    MDDir ./ utilities.rowNormAccrossM(MatrixRVec) ); % rowNormAccrossMatrix redundary
%function
eM2     = reshape(sum(construct.c0*k^3*exp(-1i*k*MR) ./ (4*pi*k*MR) .* (...
    cross((MdirR),Mp,3) + (1 ./ (k*MR*1i)+1) ),2) ,[nObs,3] );
end

