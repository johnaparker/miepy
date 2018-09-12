import miepy

S = miepy.spheroid([0,0,0], 1, 1, miepy.constant_material(2))
T = S.compute_tmatrix(2, 5, 1)
print(T)
