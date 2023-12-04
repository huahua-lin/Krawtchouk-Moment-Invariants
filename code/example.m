path = '';
dim = 128;
img = readImage(path, dim);
invariants = get3DKMI(img, dim);
