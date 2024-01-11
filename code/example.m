path = '../example_data/airplane.im';
dim = 128;
img = readImage(path, dim);
invariants1 = get3DKMI(img, dim);

%%
path = '../example_data/cup.im';
dim = 128;
img = readImage(path, dim);
invariants2 = get3DKMI(img, dim);
