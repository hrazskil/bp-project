clc
clear
nDip=2;
dip.pos = rand(nDip,3);
dip.dir = rand(nDip,3);
complAmpl = rand(nDip,1)+  1i*rand(nDip,1);

f = 10000;

nObs=5;
rObserved=rand(nObs,3);
%%
% tic
% [eZ] = fieldEvaluation.eleField(rObserved,dip,f,complAmpl);
% toc
%

%MAGFIELD Summary of this function goes here
%   Detailed explanation goes here
%inicialization
construct   = utilities.constants.giveConstants;
nObs        = size(rObserved,1);
mF          = nan(nObs,3);

%out for-loop variables
omega       = 2*pi*f;
k           = omega/construct.c0;

    for iObs = 1:nObs

    Rvec    = repmat(rObserved(iObs,:),[size(dip.pos,1),1])-dip.pos;
    R       = utilities.rowNorm(Rvec);
    dirR    = Rvec./repmat(R,[1,3]);
% not sure about the p, but the dimensions fit and matlab computes.
    p       = repmat(complAmpl,[1,3]).*( ...
        dip.dir./repmat(utilities.rowNorm(Rvec),[1,3]));

    mF(iObs,:) = sum(construct.c0*k^3*exp(-1i*k*R)./(4*pi*k*R) .*(...
        cross((dirR),p,2) + (1./(k*R*1i)+1) ) ,1);
    end


% tic
% [eM] = fieldEvaluation.magField(rObserved,dip,f,complAmpl);
% toc


% tic
% [eM] = fieldEvaluation.magField(rObserved,dip,f,complAmpl);
% toc

% tic
% [eM] = fieldEvaluation.magField(rObserved,dip,f,complAmpl);
% toc
%

%inicialization
construct   = utilities.constants.giveConstants;
nDip = size(dip.pos,1);
nObs = size(rObserved,1);
eM2          = nan(nObs,nDip,3);     % Preallocate
MatrixRVec  = zeros(nObs,nDip,3);   % Preallocate
MR           = zeros(nObs,nDip);     % Preallocate
MdirR        = zeros(nObs,nDip,3);   % Preallocate
omega       = 2*pi*f;
k           = omega/construct.c0;


MatrixRVec  = reshape( repmat(rObserved,nDip,1) ,nObs, [],3) -...
    permute( reshape( repmat(dip.pos,nObs,1) ,nDip,[],3) ,[2 1 3]);
MR           = utilities.rowNormAccrossM(MatrixRVec);
MdirR        = MatrixRVec./MR;


MatrixDDir  = permute( reshape( repmat(dip.dir,nObs,1) ,nDip,[],3) ...
    ,[2 1 3]);
Mp           = repmat(complAmpl,[1,nObs])'.*( ...
    MatrixDDir./repmat(utilities.rowNorm(dip.dir),[1,nObs])' ); % mistake

eM2  = sum(construct.c0*k^3*exp(-1i*k*MR)./(4*pi*k*MR) .*(...
    cross((MdirR),Mp,3) + (1./(k*MR*1i)+1) ),2)

