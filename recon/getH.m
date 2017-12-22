function H = getH(sigma_type,const,N)

Hij_fname = getSaveFname('Hij',sigma_type,const,N);
H_struct = load(Hij_fname);

H = zeros(3,3,N,N,N);

H(1,1,:,:,:) = H_struct.H11;
H(2,1,:,:,:) = H_struct.H21;
H(3,1,:,:,:) = H_struct.H31;

H(1,2,:,:,:) = H_struct.H12;
H(2,2,:,:,:) = H_struct.H22;
H(3,2,:,:,:) = H_struct.H32;

H(1,3,:,:,:) = H_struct.H13;
H(2,3,:,:,:) = H_struct.H23;
H(3,3,:,:,:) = H_struct.H33;

end
