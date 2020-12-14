
Required steps to make a model

1) perl learning_matrix.pl 
with parameters :
--mode train or test (build a labeled (train) or unlabeled (test) matrix depending on the use case.
--network input a gold standard e.g. Golden_standard_iGRN.txt only required with train mode
—-expr eg COE_iGRN.txt
—-matrix eg CB_motif_iGRN.txt
—-cluster eg CB_cluster_iGRN.txt
—-DH eg CB_DH_iGRN.txt
—-CNS eg CB_CMM_iGRN.txt
—-ChIP eg ChIP_scores_rank_iGRN.txt
—-GENIE3 eg GENIE3_iGRN.txt
—KNN eg CB_KNN_iGRN.txt

2)learning.R
Notebook to train a model based on a learning_matrix.txt or predict based on a generated test_matrix.txt