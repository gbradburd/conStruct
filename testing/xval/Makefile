.PHONY : test clean

test : clean
	Rscript test.xval.R &> test.log
	Rscript test.xval2.R &> test2.log

clean :
	rm  -f test1.xvals.Robj xval_test1.xval.data.partitions.Robj xval_test1.xval.results.Robj xval_test1_nsp_xval_results.txt xval_test1_sp_xval_results.txt test.log test2.log xval_test2.xval.data.partitions.Robj xval_test2.xval.results.Robj xval_test2_nsp_xval_results.txt xval_test2_sp_xval_results.txt xvals2.Robj .Rapp.history .RData