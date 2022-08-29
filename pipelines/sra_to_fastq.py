from argparse import ArgumentParser
import os
import sys

def _parse_cmdl(cmdl):
	parser = ArgumentParser(description="Automatic GEO SRA data downloader")
	
	parser.add_argument(
			"-b", "--bamfolder", dest="bam_folder",
			default=safe_echo("SRABAM"),
			help="Optional: Specify a location to store bam files "
			"[Default: $SRABAM:" + safe_echo("SRABAM") + "]")
	
	parser.add_argument(
			"-s", "--srafolder", dest="sra_folder", default=safe_echo("SRARAW"),
			help="Optional: Specify a location to store sra files "
			"[Default: $SRARAW:" + safe_echo("SRARAW") + "]")
	
	parser.add_argument(
			"--picard", dest="picard_path", default=safe_echo("PICARD"),
			help="Specify a path to the picard jar, if you want to convert "
			"fastq to bam [Default: $PICARD:" + safe_echo("PICARD") + "]")
	
	parser.add_argument(
			"-r", "--srr", required=True, nargs="+",
			help="SRR accession")

parser.add_argument(
			"-o", "--bam", required=True, nargs="+",
			help="BAM outfile")
	

	return parser.parse_args(cmdl)

def safe_echo(var):
	""" Returns an environment variable if it exists, or an empty string if not"""
	return os.getenv(var, "")


def main(cmdl):
	
	args = _parse_cmdl(cmdl)
	if len(args.srr) != len(args.bam):
		raise Exception("Numbers of inputs and outputs must match.")
	for i in range(len(args.srr)):
		infile = args.srr[i]
		outfile = args.bam[i]
		cmd = "fastq-dump --split-3 -O " + \
		os.path.realpath(args.sra_folder) + " " + \
		os.path.join(args.sra_folder, args.srr + ".sra")

		print(cmd)

if __name__ == "__main__":
	main(sys.argv[1:])
