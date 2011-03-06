file shoedata.csv
vars x x2 year id firm shoe y ife ffe sfe yfe
model y ~ x + x2 +G(firm)+G(id)+G(shoe)+G(year)
dummy year
#tol 1e-4
# nofe
merge

