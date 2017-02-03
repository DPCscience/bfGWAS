.DELETE_ON_ERROR:

all: /net/fantasia/home/yjingj/GIT/SFBA_TEST/pre_em.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.0.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.0.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.0.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.0.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param0.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/R0.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.1.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.1.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.1.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.1.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param1.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/R1.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.2.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.2.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.2.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.2.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param2.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/R2.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.3.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.3.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.3.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.3.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param3.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/R3.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.4.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.4.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.4.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.4.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param4.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/R4.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.5.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.5.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.5.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.5.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param5.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/R5.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/pre_em.OK: 
	rm -f -r /net/fantasia/home/yjingj/GIT/SFBA_TEST/output /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT
	mkdir -p /net/fantasia/home/yjingj/GIT/SFBA_TEST/output /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT
	cp -f  /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current
	> /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/EM_result.txt
	> /net/fantasia/home/yjingj/GIT/SFBA_TEST/Rout.txt
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/pre_em.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.0.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/pre_em.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /CFH_REGION_1KG.vcf.gz -vcfp  -a /Anno_CFH_REGION_1KG.gz -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o CFH_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.output.txt 
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.0.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.0.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/pre_em.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /CFI_REGION_1KG.vcf.gz -vcfp  -a /Anno_CFI_REGION_1KG.gz -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o CFI_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.output.txt 
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.0.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.0.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/pre_em.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /C2_REGION_1KG.vcf.gz -vcfp  -a /Anno_C2_REGION_1KG.gz -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o C2_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.output.txt 
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.0.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.0.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/pre_em.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /C3_REGION_1KG.vcf.gz -vcfp  -a /Anno_C3_REGION_1KG.gz -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o C3_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.output.txt 
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.0.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param0.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.0.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.0.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.0.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.0.OK 
	cat `ls -d -1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/output/** | grep paramtemp | sort` > /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/paramtemp0.txt
	cat `ls -d -1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/output/** | grep hyptemp | sort` > /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/hyptemp0.txt
	cat `ls -d -1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/output/** | grep log | sort` > /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/log0.txt
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param0.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/R0.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param0.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/pre_em.OK
	Rscript --no-save --no-restore --verbose ./bin/Mstep.r /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/hyptemp0.txt 0 1e-6 0.1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/EM_result.txt /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current >> /net/fantasia/home/yjingj/GIT/SFBA_TEST/Rout.txt
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/R0.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.1.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/R0.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /CFH_REGION_1KG.vcf.gz -a /Anno_CFH_REGION_1KG.gz -vcfp  -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o CFH_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.output.txt  
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.1.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.1.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/R0.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /CFI_REGION_1KG.vcf.gz -a /Anno_CFI_REGION_1KG.gz -vcfp  -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o CFI_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.output.txt  
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.1.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.1.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/R0.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /C2_REGION_1KG.vcf.gz -a /Anno_C2_REGION_1KG.gz -vcfp  -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o C2_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.output.txt  
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.1.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.1.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/R0.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /C3_REGION_1KG.vcf.gz -a /Anno_C3_REGION_1KG.gz -vcfp  -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o C3_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.output.txt  
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.1.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param1.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.1.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.1.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.1.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.1.OK 
	cat `ls -d -1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/output/** | grep paramtemp | sort ` > /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/paramtemp1.txt
	cat `ls -d -1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/output/** | grep hyptemp | sort` > /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/hyptemp1.txt
	cat `ls -d -1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/output/** | grep log | sort` > /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/log1.txt
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param1.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/R1.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param1.OK
	Rscript --no-save --no-restore --verbose ./bin/Mstep.r /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/hyptemp1.txt 1 1e-6 0.1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/EM_result.txt /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current >> /net/fantasia/home/yjingj/GIT/SFBA_TEST/Rout.txt
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/R1.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.2.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/R1.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /CFH_REGION_1KG.vcf.gz -a /Anno_CFH_REGION_1KG.gz -vcfp  -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o CFH_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.output.txt  
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.2.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.2.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/R1.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /CFI_REGION_1KG.vcf.gz -a /Anno_CFI_REGION_1KG.gz -vcfp  -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o CFI_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.output.txt  
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.2.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.2.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/R1.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /C2_REGION_1KG.vcf.gz -a /Anno_C2_REGION_1KG.gz -vcfp  -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o C2_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.output.txt  
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.2.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.2.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/R1.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /C3_REGION_1KG.vcf.gz -a /Anno_C3_REGION_1KG.gz -vcfp  -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o C3_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.output.txt  
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.2.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param2.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.2.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.2.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.2.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.2.OK 
	cat `ls -d -1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/output/** | grep paramtemp | sort ` > /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/paramtemp2.txt
	cat `ls -d -1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/output/** | grep hyptemp | sort` > /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/hyptemp2.txt
	cat `ls -d -1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/output/** | grep log | sort` > /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/log2.txt
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param2.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/R2.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param2.OK
	Rscript --no-save --no-restore --verbose ./bin/Mstep.r /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/hyptemp2.txt 2 1e-6 0.1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/EM_result.txt /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current >> /net/fantasia/home/yjingj/GIT/SFBA_TEST/Rout.txt
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/R2.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.3.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/R2.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /CFH_REGION_1KG.vcf.gz -a /Anno_CFH_REGION_1KG.gz -vcfp  -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o CFH_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.output.txt  
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.3.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.3.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/R2.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /CFI_REGION_1KG.vcf.gz -a /Anno_CFI_REGION_1KG.gz -vcfp  -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o CFI_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.output.txt  
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.3.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.3.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/R2.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /C2_REGION_1KG.vcf.gz -a /Anno_C2_REGION_1KG.gz -vcfp  -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o C2_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.output.txt  
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.3.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.3.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/R2.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /C3_REGION_1KG.vcf.gz -a /Anno_C3_REGION_1KG.gz -vcfp  -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o C3_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.output.txt  
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.3.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param3.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.3.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.3.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.3.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.3.OK 
	cat `ls -d -1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/output/** | grep paramtemp | sort ` > /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/paramtemp3.txt
	cat `ls -d -1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/output/** | grep hyptemp | sort` > /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/hyptemp3.txt
	cat `ls -d -1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/output/** | grep log | sort` > /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/log3.txt
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param3.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/R3.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param3.OK
	Rscript --no-save --no-restore --verbose ./bin/Mstep.r /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/hyptemp3.txt 3 1e-6 0.1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/EM_result.txt /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current >> /net/fantasia/home/yjingj/GIT/SFBA_TEST/Rout.txt
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/R3.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.4.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/R3.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /CFH_REGION_1KG.vcf.gz -a /Anno_CFH_REGION_1KG.gz -vcfp  -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o CFH_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.output.txt  
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.4.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.4.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/R3.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /CFI_REGION_1KG.vcf.gz -a /Anno_CFI_REGION_1KG.gz -vcfp  -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o CFI_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.output.txt  
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.4.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.4.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/R3.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /C2_REGION_1KG.vcf.gz -a /Anno_C2_REGION_1KG.gz -vcfp  -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o C2_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.output.txt  
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.4.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.4.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/R3.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /C3_REGION_1KG.vcf.gz -a /Anno_C3_REGION_1KG.gz -vcfp  -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o C3_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.output.txt  
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.4.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param4.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.4.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.4.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.4.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.4.OK 
	cat `ls -d -1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/output/** | grep paramtemp | sort ` > /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/paramtemp4.txt
	cat `ls -d -1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/output/** | grep hyptemp | sort` > /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/hyptemp4.txt
	cat `ls -d -1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/output/** | grep log | sort` > /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/log4.txt
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param4.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/R4.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param4.OK
	Rscript --no-save --no-restore --verbose ./bin/Mstep.r /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/hyptemp4.txt 4 1e-6 0.1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/EM_result.txt /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current >> /net/fantasia/home/yjingj/GIT/SFBA_TEST/Rout.txt
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/R4.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.5.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/R4.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /CFH_REGION_1KG.vcf.gz -a /Anno_CFH_REGION_1KG.gz -vcfp  -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o CFH_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.output.txt  
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.5.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.5.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/R4.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /CFI_REGION_1KG.vcf.gz -a /Anno_CFI_REGION_1KG.gz -vcfp  -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o CFI_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.output.txt  
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.5.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.5.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/R4.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /C2_REGION_1KG.vcf.gz -a /Anno_C2_REGION_1KG.gz -vcfp  -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o C2_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.output.txt  
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.5.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.5.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/R4.OK
	srun --exclude= --partition=nomosix --mem-per-cpu=3000 --time=24:00:00 --nice=100 --error=/net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/%N.%j.err -J  -D /net/fantasia/home/yjingj/GIT/SFBA_TEST /bin/Estep_mcmc -vcf /C3_REGION_1KG.vcf.gz -a /Anno_C3_REGION_1KG.gz -vcfp  -fcode  -hfile /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current -GTfield GT -maf 0.005 -bslmm -rmin 1 -rmax 1 -smin 0 -smax 5 -win 100 -n 1 -o C3_REGION_1KG -w 50000 -s 50000 -comp 0 -saveSNP 0 -initype 3 -rv 1 > /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.output.txt  
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.5.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param5.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFH_REGION_1KG.5.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/CFI_REGION_1KG.5.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C2_REGION_1KG.5.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/C3_REGION_1KG.5.OK 
	cat `ls -d -1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/output/** | grep paramtemp | sort ` > /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/paramtemp5.txt
	cat `ls -d -1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/output/** | grep hyptemp | sort` > /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/hyptemp5.txt
	cat `ls -d -1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/output/** | grep log | sort` > /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/log5.txt
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param5.OK

/net/fantasia/home/yjingj/GIT/SFBA_TEST/R5.OK: /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/cp_param5.OK
	Rscript --no-save --no-restore --verbose ./bin/Mstep.r /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/hyptemp5.txt 5 1e-6 0.1 /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/EM_result.txt /net/fantasia/home/yjingj/GIT/SFBA_TEST/hypval.current >> /net/fantasia/home/yjingj/GIT/SFBA_TEST/Rout.txt
	touch /net/fantasia/home/yjingj/GIT/SFBA_TEST/R5.OK

clean_err: 
	-rm -rf /net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/*.err
clean: 
	-rm -rf /net/fantasia/home/yjingj/GIT/SFBA_TEST/*.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/Eoutput/*.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/OUT/*.OK /net/fantasia/home/yjingj/GIT/SFBA_TEST/slurm_err/*.err
