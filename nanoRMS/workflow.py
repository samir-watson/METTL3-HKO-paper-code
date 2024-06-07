from gwf import Workflow
from gwf import AnonymousTarget
import itertools
import os.path
import csv

gwf = Workflow(defaults={'cores':8,'memory':'16gb','account':'pkadbioinfproject','walltime':'04:00:00'})
transcripts="data/genome/gencode.v33.transcripts.fa"
transcripts_shorthead = "data/genome/gencode.v33.transcripts.shorthead.fa"

def nanopolish_index(path):
	description = "nanopolish_index"
	inputs = gwf.glob('{path}/fastq_runid_*'.format(path=path))
	outputs = [path+"/run.fastq.gz",path+"/run.fastq.gz.index"]
	options = {}
	spec = """cat {path}/*.fastq.gz > {path}/run.fastq.gz && nanopolish index -s {path}/sequencing_summary.txt -d {path}/workspace/ {path}/run.fastq.gz """.format(path=path)
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def mapping(path):
	description = "mapping"
	inputs = [path+"/run.fastq.gz",transcripts]
	outputs = [path+"/run.aligned.bam"]
	options = {}
	spec = """minimap2 -ax map-ont -t 16 -L {transcripts} {path}/run.fastq.gz | samtools view -bh -F 2324 -q 10 | samtools sort -@16 -O bam > {path}/run.aligned.bam """.format(transcripts=transcripts,path=path)
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def samtools_index(path):
	inputs = [path+"/run.aligned.bam"]
	outputs = [path+"/run.aligned.bam.bai"]
	options = {}
	spec = """samtools index -@ 16 {path}/run.aligned.bam """.format(path=path)
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def eventalign(path):
	inputs = [path+"/run.aligned.bam",path+"/run.fastq.gz",transcripts]
	outputs = [path+"/run.transcriptome.aligned.eventalign.tsv"]
	options = {'cores':30,'walltime':'2-00:00:00'}
	spec = """nanopolish eventalign -t 30 --reads {path}/run.fastq.gz --bam {path}/run.aligned.bam --genome {transcripts} --print-read-names --scale-events --samples > {path}/run.transcriptome.aligned.eventalign.tsv""".format(transcripts=transcripts,path=path)
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def nanorms_picard():
	inputs = ["scripts/nanoRMS/epinano_RMS/picard.jar",transcripts]
	outputs = ["data/genome/gencode.v33.transcripts.fa.dict"]
	options = {}
	spec = """java -jar {picard} CreateSequenceDictionary -R {transcripts} -O {output}""".format(transcripts=transcripts,picard=inputs[0],output=outputs[0])
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def nanorms_epinano_rms(path):
	basedir="~/pkadbioinfproject/testguppy"
	inputs = [path+"/run.aligned.bam","scripts/nanoRMS/epinano_RMS","data/genome/gencode.v33.transcripts.fa.dict"]
	outputs = [path+"/run.aligned.per.site.baseFreq.csv"]
	options = {}
	spec = """ cd {path} && python {basedir}/scripts/nanoRMS/epinano_RMS/epinano_rms.py -R {basedir}/data/genome/gencode.v33.transcripts.fa -b run.aligned.bam -s {basedir}/scripts/nanoRMS/epinano_RMS/sam2tsv.jar""".format(basedir=basedir,path=path)
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)	
# 4 hours is def not enough for very large files

def nanorms_get_features(path):
	inputs = gwf.glob('{path}/workspace/FA*.fast5'.format(path=path)) + [transcripts_shorthead]
	
	
	outputs = [f(x) for x in inputs[:-1] for f in (lambda x:x+".bam",lambda x:x+".bam.json",lambda x:x+".bam.bai")]

	options = {"cores":16,"memory":"128gb","walltime":"6:00:00"}
	spec = """
	python scripts/nanoRMS/per_read/get_features.py --rna -f {transcripts} -t 16 -i {path}/workspace/FA*.fast5
	""".format(transcripts=transcripts_shorthead,path=path)
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def nanorms_per_read_mean(path):
	inputs = [path+"/run.transcriptome.aligned.eventalign.tsv"]
	outputs = [path+"/run.transcriptome.aligned.eventalign.tsv_processed_perpos_mean.tsv"]
	options = {"memory":"256gb","walltime":"8:00:00"}
	spec = """
	python3 scripts/nanoRMS/visualization_per_read/Dper_read_mean.py {path}/run.transcriptome.aligned.eventalign.tsv

	""".format(path=path)
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def nanorms_per_read_mean_fixcontig(path):
	inputs = [path+"/run.transcriptome.aligned.eventalign.tsv_processed_perpos_mean.tsv"]
	outputs = [path+"/run.transcriptome.aligned.eventalign.tsv_processed_perpos_mean_fixedcontig.tsv"]
	options = {"cores":2,"memory":"4gb","walltime":"3:00:00"} # probably dosnt need more than a gb of memory really
	spec = """
	{command} {path}/run.transcriptome.aligned.eventalign.tsv_processed_perpos_mean.tsv > {path}/run.transcriptome.aligned.eventalign.tsv_processed_perpos_mean_fixedcontig.tsv
	""".format(path=path,command=r"""awk -F '\t' 'BEGIN {OFS=FS} {split($1,a,"|"); $1=a[1]; print}'""") #kinda weird but circumvents escape hell in the awk command, i hope
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def nanorms_plots(path): #dummy function spec is from fixcontig
	inputs = [path+""]
	outputs = [path+"/run.transcriptome.aligned.eventalign.tsv_processed_perpos_mean_fixedcontig.tsv"]
	options = {"memory":"64gb","walltime":"2:00:00"}
	spec = """
	{command} {path}/run.transcriptome.aligned.eventalign.tsv_processed_perpos_mean.tsv > {path}/run.transcriptome.aligned.eventalign.tsv_processed_perpos_mean_fixedcontig.tsv
	""".format(path=path,command=r"""awk -F '\t' 'BEGIN {OFS=FS} {split($1,a,"|"); $1=a[1]; print}'""") #kinda weird but circumvents escape hell in the awk command, i hope
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def nanoqc(path):
	inputs = [path+"/run.aligned.tmp_splitted_base_freq"]
	outputs = []
	options = {}
	spec = """
	per_read/get_freq.py -f {transcripts} -b $f.bed -o $f.bed.tsv.gz -1 per_read/guppy3.0.3.hac/*WT30C/workspace/*.fast5.bam -2 per_read/guppy3.0.3.hac/*WT45C/workspace/*.fast5.bam
	""".format(transcripts=transcripts,)
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def window(contig,position,label): # memory needs a lot on first run when creating parquet, and half on subsequent runs. maybe figure out logic for this?
	description = "kmer windows"
	inputs = ["data/METTL3/{label}_HQcalls_RMS/run.transcriptome.aligned.eventalign.tsv_processed_perpos_mean_fixedcontig.tsv".format(label=label)]
	outputs = ["steps/windows/{label}_windows/{label}_{contig}_{position}.window_file.tsv.tmp".format(label=label,contig=contig,position=position)]
	options = {'cores':2,'memory':'110gb','account':'pkadbioinfproject','walltime':'02:00:00'}
	spec = """python scripts/windowscript.py {contig} {position} {input} {label} """.format(contig=contig,position=position, input=inputs[0], label=label)
	
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def window15mer(contig,position,label): 
	description = "kmer windows perread"
	inputs = ["steps/windows/{label}_windows/{label}_{contig}_{position}.window_file.tsv.tmp".format(label=label,contig=contig,position=position)]
	outputs = ["steps/windows15mer/{label}_windows/{contig}_{position}_{label}_15mer.perread.tsv".format(label=label,contig=contig,position=position)]
	options = {'cores':1,'memory':'5gb','account':'pkadbioinfproject','walltime':'00:20:00'}
	spec = """Rscript --vanilla scripts/nanoRMS/visualization_per_read/nanopolish_export_each_position_singleinput.R {input} {output}/""".format(input=inputs[0],output=os.path.dirname(outputs[0]))
	
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def clustering(contig,position): #Hardcoded to WT_2 and A56_2
	description = "kmeans clustering of ENST position"
	inputs = ["steps/windows15mer/WT_2_windows/{contig}_{position}_WT_2_15mer.perread.tsv".format(contig=contig,position=position),"steps/windows15mer/A56_2_windows/{contig}_{position}_A56_2_15mer.perread.tsv".format(contig=contig,position=position)]
	
	outputs = [f("results/clustering/{contig}/{position}/".format(contig=contig,position=position)) for f in (lambda x:x+"kmeans_colorBY_SAMPLE.pdf",lambda x:x+"kmeans_colorBY_predictedModStatus.pdf",lambda x:x+"kmeans.txt",lambda x:x+"knn.txt")]
	
	options = {'cores':1,'memory':'5gb','account':'pkadbioinfproject','walltime':'00:20:00'}
	spec = """R --vanilla < scripts/nanoRMS/visualization_per_read/read_clustering_fixed.R --args {input} kmeans && R --vanilla < scripts/nanoRMS/visualization_per_read/read_clustering_fixed.R --args {input} knn""".format(input=contig+"_"+position)
	
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def tombodiffmod(positions):
	description = "kmeans clustering of ENST position"
	WT = gwf.glob("data/METTL3/WT_2_HQcalls_RMS/workspace/*.fast5.bam")
	KO = gwf.glob("data/METTL3/A56_2_HQcalls_RMS/workspace/*.fast5.bam")
	inputs = WT + KO  + ["/home/peterkad/pkadbioinfproject/testguppy/steps/nanocomporeresult_chr_nodrach_ENST_withstrand.tsv",transcripts_shorthead]
	outputs = ["results/tombo_diffmod_{positions}".format(positions=os.path.basename(positions))]
	
	options = {'cores':4,'memory':'100gb','account':'pkadbioinfproject','walltime':'2-12:00:00'}

	spec = """python scripts/nanoRMS/per_read/get_freq_modified.py -f {transcripts_shorthead} -b {positions} -o {outputfile} -1 {wt_input} -2 {ko_input}""".format(transcripts_shorthead=transcripts_shorthead,positions=positions,outputfile=outputs[0],wt_input="data/METTL3/WT_2_HQcalls_RMS/workspace/*.fast5.bam",ko_input="data/METTL3/A56_2_HQcalls_RMS/workspace/*.fast5.bam")
	
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def nanopolish_plots(contig,position):
	description = "plots current intensity distribution, mean current intensity centered in window and PCA plots"
	
	inputs = ["steps/windows/{label}_windows/{label}_{contig}_{position}.window_file.tsv.tmp".format(label=label,contig=contig,position=position) for label in ["WT_2","A56_2"]]
	outputs = ["plots/nanopolish/{contig}_{position}_{plot}".format(contig=contig,position=position,plot=plot) for plot in ["mean_line_plot.naomit.pdf","pca.pdf","density_plot.pdf","perread_lineplot.pdf"]]
	
	options = {'cores':1,'memory':'16gb','account':'pkadbioinfproject','walltime':'1:00:00'}
	spec = """ Rscript --vanilla scripts/nanoRMS/visualization_per_read/nanopolish_density_plot.R {wt} {A56} && Rscript --vanilla scripts/nanoRMS/visualization_per_read/nanopolish_meanlineplot.R {wt} {A56} && Rscript --vanilla scripts/nanoRMS/visualization_per_read/nanopolish_pca.R {wt} {A56} && Rscript --vanilla scripts/nanoRMS/visualization_per_read/nanopolish_perreadlineplot.R  {wt} {A56} && mv ENST* plots/nanopolish/ """.format(wt=inputs[0],A56=inputs[1])
	
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


### this is stupid, there has to be a gooder way.
def get_name_nanopolish_index(idx,target):
	filename = os.path.basename(os.path.dirname(target.inputs[0]))
	return '{function}_{filename}'.format(filename=filename,function="nanopolish_index")
def get_name_mapping(idx,target):
	filename = os.path.basename(os.path.dirname(target.inputs[0]))
	return '{function}_{filename}'.format(filename=filename,function="mapping")
def get_name_samtools_index(idx,target):
	filename = os.path.basename(os.path.dirname(target.inputs[0]))
	return '{function}_{filename}'.format(filename=filename,function="samtools_index")
def get_name_eventalign(idx,target):
	filename = os.path.basename(os.path.dirname(target.inputs[0]))
	return '{function}_{filename}'.format(filename=filename,function="eventalign")
def get_name_nanorms_epinano_rms(idx,target):
	filename = os.path.basename(os.path.dirname(target.inputs[0]))
	return '{function}_{filename}'.format(filename=filename,function="nanorms_epinano_rms")
def get_name_nanorms_get_features(idx,target):
	filename = os.path.basename(os.path.split(os.path.dirname(target.inputs[0]))[0])
	return '{function}_{filename}'.format(filename=filename,function="nanorms_get_features")	
def get_name_nanorms_per_read_mean(idx,target):
	filename = os.path.basename(os.path.dirname(target.inputs[0]))
	return '{function}_{filename}'.format(filename=filename,function="nanorms_per_read_mean")
def get_name_nanorms_per_read_mean_fixcontig(idx,target):
	filename = os.path.basename(os.path.dirname(target.inputs[0]))
	return '{function}_{filename}'.format(filename=filename,function="nanorms_per_read_mean_fixcontig")
def get_name_nanorms_window15mer(idx,target):
	filename = os.path.basename(os.path.dirname(target.inputs[0]))
	return '{function}_{filename}'.format(filename=filename,function="nanorms_window15mer")			


## mapping targets and executing the workflow

targets = gwf.glob('data/METTL3/*_?_HQcalls_RMS')
gwf.map(nanopolish_index, targets,name=get_name_nanopolish_index)
gwf.map(mapping, targets,name=get_name_mapping)
gwf.map(samtools_index, targets,name=get_name_samtools_index)
gwf.map(eventalign, targets,name=get_name_eventalign)
gwf.target_from_template(name="reference_dict",template=nanorms_picard())
gwf.map(nanorms_epinano_rms,targets,name=get_name_nanorms_epinano_rms)
gwf.map(nanorms_get_features,targets,name=get_name_nanorms_get_features)
gwf.map(nanorms_per_read_mean,targets,name=get_name_nanorms_per_read_mean)
gwf.map(nanorms_per_read_mean_fixcontig,targets,name=get_name_nanorms_per_read_mean_fixcontig)

#replicate labels 
#,"A56_2" A56_2 is absolutly terrible and should be discarded.
#labels = ["WT_1","WT_2","WT_3","A56_1","A56_2","A56_3","A56_1","A56_3"]
labels = ["WT_2","A56_2"]
if False: # this is roughly 1700 targets per replicate. This kills gwf if run every time.
	with open('/home/peterkad/pkadbioinfproject/testguppy/steps/nanocomporeresult_chr_nodrach_ENST.tsv', 'r') as f:
		lines = f.read().splitlines()[1:] # discard the first line
		contig = [line.split('\t')[0] for line in lines]
		pos = [line.split('\t')[1] for line in lines]


	for contig,pos in zip(contig,pos):
		for label in labels:
			gwf.target_from_template('window_{}_{}_{}'.format(label,contig, pos), window(contig,pos,label))

def windowOfInterest_runner(input): 
	interests = [(row[0],row[1]) for row in csv.reader(open(input,'r'), delimiter='\t')][1:]
	for windowI in ((touple + (label,)) for touple in interests for label in labels):
		gwf.target_from_template('window_15mer_{}_{}_{}'.format(*windowI), window15mer(*windowI))
	return None

#windowOfInterest_runner('/home/peterkad/pkadbioinfproject/testguppy/steps/nanocomporeresult_sorted2.tsv')
#windowOfInterest_runner('/home/peterkad/pkadbioinfproject/testguppy/steps/nanocomporeresult_sorted_DRACH_unique.tsv')

def clusterOfInterest(input):
	interests = [(row[0],row[1]) for row in csv.reader(open(input,'r'), delimiter='\t')][1:]
	for clusterI in interests:
		gwf.target_from_template('clustering_{}_{}'.format(*clusterI), clustering(*clusterI))
	return None

#clusterOfInterest('/home/peterkad/pkadbioinfproject/testguppy/steps/nanocomporeresult_sorted2.tsv')
#clusterOfInterest('/home/peterkad/pkadbioinfproject/testguppy/steps/nanocomporeresult_sorted_DRACH_unique.tsv')

#gwf.target_from_template('tombodiffmod_WT_2VSA56_2',tombodiffmod("steps/nanocomporeresult_chr_nodrach_ENST_withstrand.tsv"))

gwf.target_from_template('tombodiffmod_WT_2VSA56_2_MIFCYBMTRNR',tombodiffmod("steps/nanocomporeresult_chr_nodrach_ENST_withstrandCYBMIFMTRNR.tsv"))
gwf.target_from_template('tombodiffmod_WT_2VSA56_2_top100',tombodiffmod("steps/nanocomporeresult_chr_nodrach_ENST_withstrand_top100.tsv"))

for contig,position in [("ENST00000387347.2","1318"),("ENST00000361789.2","213"),("ENST00000215754.8", "281"),("ENST00000260379.11", "130"),("ENST00000362079.2", "69"),("ENST00000270625.7", "504")]:
	gwf.target_from_template(f"nanopolish_plots_{contig}_{position}",nanopolish_plots(contig,position))
	gwf.target_from_template(f'clustering_{contig}_{position}', clustering(contig,position))
	gwf.target_from_template('window_15mer_{}_{}_{}'.format(contig,position,labels[0]), window15mer(contig,position,labels[0]))
	gwf.target_from_template('window_15mer_{}_{}_{}'.format(contig,position,labels[1]), window15mer(contig,position,labels[1]))

#gwf.target_from_template("clustering_ENST00000215754.8_515", clustering("ENST00000215754.8","515"))
#gwf.target_from_template('window2fullframe', window("ENST00000007516.8","490","A56_2"))
#gwf.target_from_template('window2pdhdf', window("ENST00000007516.8","490","A56_3"))