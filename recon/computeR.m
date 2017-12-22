function R = computeR(aa,bb,cc,dd)

R = [aa^2 + bb^2 - cc^2 - dd^2, 2*(bb*cc-aa*dd), 2*(bb*dd+aa*cc); ...
    2*(bb*cc + aa*dd), aa^2 - bb^2 + cc^2 - dd^2, 2*(cc*dd-aa*bb); ...
    2*(bb*dd-aa*cc), 2*(cc*dd+aa*bb), aa^2 - bb^2 - cc^2 + dd^2];

end