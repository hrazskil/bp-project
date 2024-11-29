function [powerFull] = powerFull(eF,mF)
%POWERFULL computes power of electromagnetic field in submited field points
powerFull = real(cross( eF, conj(mF),2))/2;
end

