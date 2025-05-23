function [encoded, trellis] = FEC_Coding(data)
dataT = reshape(data,[],1);
% trellis = poly2trellis(4,{'1+x+x^2+x^3','1+x+x^3'});
trellis = poly2trellis(4,[17, 13]);
codedDataT = convenc(dataT,trellis);
% codedData = reshape(codedDataT,[],1);

encoded = codedDataT;


end
