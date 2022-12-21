rm -rf bin
mkdir bin
cd src
g++-8 -std=c++17 preprocess_main.cpp preprocess.cpp -o preprocess
cp preprocess ../bin
g++-8 -std=c++17 sg_lasso_main.cpp sg_lasso.cpp -o sg_lasso -Iinclude -lopenblas -llapack
cp sg_lasso ../bin
g++-8 -std=c++17 sg_lasso_leastr_main.cpp sg_lasso_leastr.cpp -o sg_lasso_leastr -Iinclude -lopenblas -llapack
cp sg_lasso_leastr ../bin
g++-8 -std=c++17 overlapping_sg_lasso_leastr_main.cpp overlapping_sg_lasso_leastr.cpp -o overlapping_sg_lasso_leastr -Iinclude -lopenblas -llapack
cp overlapping_sg_lasso_leastr ../bin
g++-8 -std=c++17 overlapping_sg_lasso_logisticr_main.cpp overlapping_sg_lasso_logisticr.cpp -o overlapping_sg_lasso_logisticr -Iinclude -lopenblas -llapack
cp overlapping_sg_lasso_logisticr ../bin
cd ..
