rm -rf angiosperm*input
rm angiosperm*feature_weights*.xml
cd sample_files
../mlpack-3.2.2/build/bin/preprocess angiosperm_20spec_pred.txt angiosperm_100_sample_alns.txt angiosperm_input
mv angiosperm_input ..
../mlpack-3.2.2/build/bin/preprocess angiosperm_20spec_pred.txt angiosperm_100_sample_alns_overlapping.txt angiosperm_ol_input
mv angiosperm_ol_input ..
cd ..
mlpack-3.2.2/build/bin/mlpack_sg_lasso_leastr -v -f angiosperm_input/feature_angiosperm_input.txt -z 0.1 -y 0.5 -n angiosperm_input/group_indices_angiosperm_input.txt -r angiosperm_input/response_angiosperm_input.txt -w angiosperm_out_feature_weights.xml
mlpack-3.2.2/build/bin/mlpack_overlapping_sg_lasso_leastr -v -f angiosperm_input/feature_angiosperm_input.txt -z 0.1 -y 0.5 -n angiosperm_input/group_indices_angiosperm_input.txt -g angiosperm_input/field_angiosperm_input.txt -r angiosperm_input/response_angiosperm_input.txt -w angiosperm_out_feature_weights_ol_fake.xml
mlpack-3.2.2/build/bin/mlpack_overlapping_sg_lasso_leastr -v -f angiosperm_ol_input/feature_angiosperm_ol_input.txt -z 0.1 -y 0.5 -n angiosperm_ol_input/group_indices_angiosperm_ol_input.txt -g angiosperm_ol_input/field_angiosperm_ol_input.txt -r angiosperm_ol_input/response_angiosperm_ol_input.txt -w angiosperm_ol_out_feature_weights_ol_real.xml
