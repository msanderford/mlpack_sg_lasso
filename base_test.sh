rm -rf angiosperm*input
rm angiosperm*feature_weights*.xml
cd sample_files
../bin/preprocess angiosperm_20spec_pred.txt angiosperm_100_sample_alns.txt angiosperm_input
mv angiosperm_input ..
../bin/preprocess angiosperm_20spec_pred.txt angiosperm_100_sample_alns_overlapping.txt angiosperm_ol_input
mv angiosperm_ol_input ..
cd ..
bin/sg_lasso -f angiosperm_input/feature_angiosperm_input.txt -z 0.1 -y 0.5 -n angiosperm_input/group_indices_angiosperm_input.txt -r angiosperm_input/response_angiosperm_input.txt -w angiosperm_out_feature_weights_logistic
bin/sg_lasso_leastr -f angiosperm_input/feature_angiosperm_input.txt -z 0.1 -y 0.5 -n angiosperm_input/group_indices_angiosperm_input.txt -r angiosperm_input/response_angiosperm_input.txt -w angiosperm_out_feature_weights
bin/overlapping_sg_lasso_leastr -f angiosperm_input/feature_angiosperm_input.txt -z 0.1 -y 0.5 -n angiosperm_input/group_indices_angiosperm_input.txt -g angiosperm_input/field_angiosperm_input.txt -r angiosperm_input/response_angiosperm_input.txt -w angiosperm_out_feature_weights_ol_fake
bin/overlapping_sg_lasso_leastr -f angiosperm_ol_input/feature_angiosperm_ol_input.txt -z 0.1 -y 0.5 -n angiosperm_ol_input/group_indices_angiosperm_ol_input.txt -g angiosperm_ol_input/field_angiosperm_ol_input.txt -r angiosperm_ol_input/response_angiosperm_ol_input.txt -w angiosperm_ol_out_feature_weights_ol_real
bin/overlapping_sg_lasso_logisticr -f angiosperm_input/feature_angiosperm_input.txt -z 0.1 -y 0.5 -n angiosperm_input/group_indices_angiosperm_input.txt -g angiosperm_input/field_angiosperm_input.txt -r angiosperm_input/response_angiosperm_input.txt -w angiosperm_out_feature_weights_ol_log_fake
bin/overlapping_sg_lasso_logisticr -f angiosperm_ol_input/feature_angiosperm_ol_input.txt -z 0.1 -y 0.5 -n angiosperm_ol_input/group_indices_angiosperm_ol_input.txt -g angiosperm_ol_input/field_angiosperm_ol_input.txt -r angiosperm_ol_input/response_angiosperm_ol_input.txt -w angiosperm_ol_out_feature_weights_ol_log_real
