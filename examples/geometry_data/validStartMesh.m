function [elements,elementgeneration,coordinates,varargout] = validStartMesh(elements,coordinates,varargin)
nE = size(elements,1);
nC = size(coordinates,1);
nB = nargin-2;
nBE = zeros(1,nB);
for j = 1:nB
    nBE(j) = size(varargin{j},1);
end
%*** Erzeugen der Flaecheninformationen
[face2nodes,element2faces,boundary2faces{1:nB}]=provideFaceData(elements,varargin{:});
nF = size(face2nodes,1);
%*** Generierung der Schwerpunkte
elementcenter = 1/4*(  coordinates(elements(:,1),:) + coordinates(elements(:,2),:) ...
                     + coordinates(elements(:,3),:) + coordinates(elements(:,4),:) );
facecenter = 1/3*(  coordinates(face2nodes(:,1),:) + coordinates(face2nodes(:,2),:) ...
                  + coordinates(face2nodes(:,3),:) );
coordinates = [coordinates;elementcenter;facecenter];
%*** Bestimmung der neuen Knotennummern
element2newNodes = nC+(1:nE)';
face2newNodes = nE+nC+(1:nF)';
%*** Erstellung der neuen Elemente
elements = [face2nodes(element2faces,1),face2newNodes(element2faces,:),repmat(element2newNodes,4,1),face2nodes(element2faces,2);...
        face2nodes(element2faces,1),face2newNodes(element2faces,:),repmat(element2newNodes,4,1),face2nodes(element2faces,3);...
        face2nodes(element2faces,2),face2newNodes(element2faces,:),repmat(element2newNodes,4,1),face2nodes(element2faces,3)];
elementgeneration = 2*ones(12*nE,1);
%*** Erstellung der Randflaechen
for j = 1:min(nargout-2,nB)
    if ~isempty(varargin{j})
        varargout{j} = [face2nodes(boundary2faces{j},1),face2newNodes(boundary2faces{j}),face2nodes(boundary2faces{j},2);...
                        face2nodes(boundary2faces{j},1),face2newNodes(boundary2faces{j}),face2nodes(boundary2faces{j},3);...
                        face2nodes(boundary2faces{j},2),face2newNodes(boundary2faces{j}),face2nodes(boundary2faces{j},3)];
    else
        varargout{j} = [];
    end
end
