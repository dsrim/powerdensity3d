function A = getMatrixAtj(Z,j)

A = [ Z.m11(j), Z.m12(j), Z.m13(j);
      Z.m21(j), Z.m22(j), Z.m23(j);
      Z.m31(j), Z.m32(j), Z.m33(j)];

end
