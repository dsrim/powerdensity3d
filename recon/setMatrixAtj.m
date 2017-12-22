function Z = setMatrixAtj(Z,A,j)

    Z.m11(j) = A(1,1);
    Z.m21(j) = A(2,1);
    Z.m31(j) = A(3,1);
    Z.m12(j) = A(1,2);
    Z.m22(j) = A(2,2);
    Z.m32(j) = A(3,2);
    Z.m13(j) = A(1,3);
    Z.m23(j) = A(2,3);
    Z.m33(j) = A(3,3);

end
