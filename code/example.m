path = '1gco.stl';
dim = 128;
img = readImage(path, dim);
invariants = get3DKMI(img, dim, 5);