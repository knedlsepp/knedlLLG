function ani = anisotropy(m,material)
%ANISOTROPY    Computes the anisotropy DPhi(m).
%   ANI = ANISOTROPY(M, MATERIAL) returns the anisotropy DPhi(M) = -2(e,x)e
%   for given magnetization m and struct MATERIAL containing the easyaxis
%   MATERIAL.easyaxis
%
%   Author: Josef Kemetmueller - 16.12.2013
e = reshape(material.easyaxis,3,[]);
ani = -2*bsxfun(@times,e',m*e);