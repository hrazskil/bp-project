clc
clear

dip.pos = rand(6,3);
dip.dir = rand(6,3);

complAmpl = rand(6,1)+  1i*rand(6,1);
f = 10000;

rObserved=rand(10,3);



%% 

% tic
% [eZ] = fieldEvaluation.eleField(rObserved,dip,f,complAmpl);
% toc
% 
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
% %inicialization
% construct   = utilities.constants.giveConstants;
% nObs        = size(rObserved,1);
% eZ          = nan(nObs,3);
% 
% %out for-loop variables
% omega       = 2*pi*frequency;
% k           = omega/construct.c0;
% 
% 
%     Rvec    = repmat(rObserved(iObs,:),[size(dip.pos,1),1])-dip.pos;
%     R       = utilities.rowNorm(Rvec);
% 
% % not sure about the p, but the dimensions fit and matlab computes.
%     p       = repmat(complAmpl,[1,3]).* ...
%         (dip.dir./repmat(utilities.rowNorm(Rvec),[1,3]));
% 
%     eZ(iObs,:) = sum(construct.Z0*construct.c0*k^3* ...
%     exp(-1i*k*R)./(4*pi*k*R) .*( cross(-(dirR),cross((dirR),p,2),2) + ...
%     ( (3*(dirR).*(dot(p,dirR,2))-p) .* (1./(k^2*R.^2)+(1i./(k*R))) ) ),1);


powerRad = fieldEvaluation.powerRadiated(f,complAmpl);

