function [arfiImage]=genArfiImg(axial, t, arfidata,normCurve)
[t2 arfidata2]=filtArfiData(axial, t, arfidata, [50 1.8e3]);
startIndex=find(t2>0,1);
arfiImage=zeros(size(arfidata2));
for i=startIndex:size(normCurve,2)
    arfidataTemp=squeeze(arfidata2(:,:,i));
    normTemp=repmat(normCurve(:,i),[1 size(arfidata2,2)]);
    arfiImage(:,:,i)=medfilt2((arfidataTemp./normTemp)*max(normTemp(end/4:(3*end/4))),[33 5]);
end

end