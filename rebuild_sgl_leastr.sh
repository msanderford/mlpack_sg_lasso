cp -R sg_lasso_src/src mlpack-3.2.2
cd mlpack-3.2.2/build
CC="gcc" cmake ../
make mlpack_sg_lasso_leastr
cd ../..
