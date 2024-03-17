import os
import subprocess

cohorts = os.listdir("cohorts")
atlas_path = "signatures/moss" #  "signatures/sinature_withoutGBM.csv"

def get_arrayTypes(idat_folder):
	cmd = "Rscript separate_arrayTypes.r cohort results"
	savedir = os.path.join("results", "cohort")
	if not os.path.exists(savedir):
		os.mkdir(savedir)
	subprocess.run(cmd, shell=True)
	df = pd.read_csv(os.path.join(savedir, "arrayTypes.csv"), index_col=0)
	arrayTypes = df.arrayTypes.unique()
	
	return df, arrayTypes
	
def process(idat_folder):
	savepath = os.path.join("results", "beta_"+idat_folder.split("/")[-1]+".csv")
	cmd = "Rscript process_array.R idat_folder savepath ref_sample.RData"
	subprocess.run(cmd, shell=True)
	
to_process = []
for cohort in cohorts:
	idat_folder = os.path.join(cohorts, cohort)
	df, arrayTypes = get_arrayTypes(idat_folder)
	if len(arrayTypes)>1:
		import shutil
		print(f"More than 1 array type in {idat_folder}.")
		for arrayType in arrayTypes:
			sub = df[df.arrayTypes==arrayType]
			tmpdir = os.path.join(cohorts, cohort+f"_{arrayType}")
			
			for f in sub.files:
				shutil.copyfile(f,  os.path.join(tmpdir, f.split("/")[-1]))
		
			to_process.append(tmp_dir)
	else:
		to_process.append(idat_folder)
		
for i in to_process:
	process(i)
	beta_file = os.path.join("results", "beta_"+idat_folder.split("/")[-1]+".csv")
	savepath = os.path.join("results", "deconvolution_"+idat_folder.split("/")[-1]+".csv")
	residual = False
	Deconvolve(atlas_path, beta_file, save_path, residual).run()

		
			
				
	
	
