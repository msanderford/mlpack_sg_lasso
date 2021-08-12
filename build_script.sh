wget https://www.mlpack.org/files/mlpack-3.2.2.tar.gz
tar -xvzpf mlpack-3.2.2.tar.gz
cp -R sg_lasso_src/src mlpack-3.2.2
cd mlpack-3.2.2
mkdir build
cd build
cmake ../
make mlpack_sg_lasso mlpack_sg_lasso_leastr mlpack_overlapping_sg_lasso_leastr
cd ../../ol_sg_lasso_preprocessing_src
g++ -std=c++11 -o preprocess preprocess_main.cpp preprocess.cpp
cd ..
cp ol_sg_lasso_preprocessing_src/preprocess mlpack-3.2.2/build/bin
